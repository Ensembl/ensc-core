/*
 * Copyright [1999-2017] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *      http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "RefineSolexaGenes.h"
#include <stdio.h>

#include <sysexits.h>
#include "EnsC.h"

#include "DBAdaptor.h"
#include "CoordSystemAdaptor.h"
#include "CachingSequenceAdaptor.h"
#include "AnalysisAdaptor.h"
#include "GeneAdaptor.h"
#include "DNAAlignFeatureAdaptor.h"
#include "BaseAdaptor.h"
#include "Vector.h"
#include "Slice.h"
#include "Species.h"
#include "StrUtil.h"
#include "DNAAlignFeature.h"
#include "DNAPepAlignFeature.h"
#include "Gene.h"
#include "Intron.h"
#include "IntronSupportingEvidence.h"
#include "Transcript.h"
#include "translate.h"
#include "Attribute.h"
#include "MetaContainer.h"

#include "libconfig.h"
#ifdef HAVE_LIBTCMALLOC
#include "gperftools/tcmalloc.h"
#endif

#include "htslib/kstring.h"
#include "htslib/sam.h"
#include "htslib/hts.h"
/*
=head1 DESCRIPTION

This module takes intron spanning dna_align_features and combines them with 
rough transcript models to build refined genes with CDS. The module produces 
transcripts representing all possible combinations of introns and exons which
are then filtered according to the specifications in the config.
The databases containing the various features to combine is defined in
Bio::EnsEMBL::Analysis::Config::Databases and the configuration for the 
module is defined in Bio::EnsEMBL::Analysis::Config::GeneBuild::RefineSolexaGenes
*/

// Used for counting recursive calls in processPaths and processTree
static int limit = 0;

static int nExonClone = 0;

#define EXIT_MEMORY 11
#define NO_GENES_WRITTEN 42

FILE *logfp;

#define RSGVERSION "0.3.2"
#define RSGGENE_KEEP 16
#define RSG_DUPLICATE 32
int dumpGenes(Vector *genes, int withSupport);

/*
typedef enum CigarBlockTypeEnum {
  CB_NONE,
  CB_MATCH,
  CB_INTRON,
  CB_DELETION
} CigarBlockType;

typedef struct CigarBlockStruct {
  CigarBlockType type;
  long  start;    // start pos on reference, 1 based I think
  long  end;      // end pos on reference, 1 based I think
} CigarBlock;
*/

int pointer_compFunc(const void *a, const void *b) {
  void *one = *((void **)a);
  void *two = *((void **)b);

  if (one == two) {
    return 0;
  } else if (one > two) {
    return 1;
  } else {
    return -1;
  }
}

CigarBlock *CigarBlock_new(CigarBlockType type, long start, long end) {
  CigarBlock *block;

  if ((block = (CigarBlock *)calloc(1,sizeof(CigarBlock))) == NULL) {
    fprintf(stderr,"Failed allocating CigarBlock\n");
    exit(EXIT_MEMORY);
  }
  block->type  = type;
  block->start = start;
  block->end   = end;

  return block;
}

CigarBlock *CigarBlock_fill(CigarBlock *block, CigarBlockType type, long start, long end) {
  block->type  = type;
  block->start = start;
  block->end   = end;

  return block;
}

void CigarBlock_copy(CigarBlock *to, CigarBlock *from) {
  memcpy(to,from,sizeof(CigarBlock));
/*
  to->type  = from->type;
  to->start = from->start;
  to->end   = from->end;
*/

  return;
}

int CigarBlock_startCompFunc(const void *a, const void *b) {
  CigarBlock *cb1 = *((CigarBlock **)a);
  CigarBlock *cb2 = *((CigarBlock **)b);

  return cb1->start - cb2->start;
}

void CigarBlock_free(CigarBlock *block) {
  free(block);
}

IntronBamConfig *IntronBamConfig_new(char *fileName, int mixedBam, int depth, Vector *groupNames) {
  IntronBamConfig *ibc;

  if ((ibc = (IntronBamConfig *)calloc(1,sizeof(IntronBamConfig))) == NULL) {
    fprintf(stderr,"Failed allocating IntronBamConfig\n");
    exit(EXIT_MEMORY);
  }
  StrUtil_copyString(&ibc->fileName, fileName, 0);
  ibc->mixedBam   = mixedBam;
  ibc->depth      = depth;
  ibc->groupNames = groupNames;

  return ibc;
}

void IntronBamConfig_free(IntronBamConfig *ibc) {
  if (ibc->fileName) free(ibc->fileName);
  
  if (ibc->groupNames) {
    Vector_setFreeFunc(ibc->groupNames, free);
    Vector_free(ibc->groupNames);
  }
  free(ibc);
}

typedef struct ModelStruct {
  Vector *features;
  double totalScore;
  char *exonUse;
} Model;

int Model_reverseScoreCompFunc(const void *a, const void *b) {
  Model *mod1 = *((Model **)a);
  Model *mod2 = *((Model **)b);

  return mod2->totalScore - mod1->totalScore;
}

int Model_firstExonStartCompFunc(const void *a, const void *b) {
  Model *mod1 = *((Model **)a);
  Model *mod2 = *((Model **)b);

  SeqFeature *sf1 = Vector_getElementAt(mod1->features, 0);
  SeqFeature *sf2 = Vector_getElementAt(mod2->features, 0);
  return SeqFeature_getStart(sf1) - SeqFeature_getStart(sf2);
}

void ModelFeatures_freeFunc(SeqFeature *sf) {
  if (sf->objectType == CLASS_EXON) {
    Exon_free(sf);
  }
}
void Model_free(Model *model) {
  if (model->exonUse) free(model->exonUse);
  
  if (model->features) {
//    Vector_setFreeFunc(model->features, ModelFeatures_freeFunc);
    Vector_free(model->features);
  }
  free(model);
}

typedef struct TranscriptExtraDataStruct {
  int intronCount;
  int nNCIntrons;
  char *exonUse;
  int depth;
  Vector *introns;
} TranscriptExtraData;

TranscriptExtraData *TranscriptExtraData_new() {
  TranscriptExtraData *ted;

  if ((ted = (TranscriptExtraData *)calloc(1,sizeof(TranscriptExtraData))) == NULL) {
    fprintf(stderr,"Failed allocating TranscriptExtraData\n");
    exit(EXIT_MEMORY);
  }

  return ted;
}

void TranscriptExtraData_free(TranscriptExtraData *ted) {
  Vector_free(ted->introns);
  free(ted->exonUse);
  free(ted);
}


ExtraExonData *ExtraExonData_new(long *coords, int nCoord) {
  ExtraExonData *eed;

  if ((eed = (ExtraExonData *)calloc(1,sizeof(ExtraExonData))) == NULL) {
    fprintf(stderr,"Failed allocating ExtraExonData\n");
    exit(EXIT_MEMORY);
  }

  eed->nCoord = nCoord;

  if ((eed->coords = (long *)calloc(nCoord,sizeof(long))) == NULL) {
    fprintf(stderr,"Failed allocating ExtraExonData coords for %d coords\n", nCoord);
    exit(EXIT_MEMORY);
  }

  memcpy(eed->coords, coords, nCoord*sizeof(long));

  return eed;
}

int ExtraExonData_startCompFunc(const void *a, const void *b) {
  ExtraExonData *eed1 = *((ExtraExonData **)a);
  ExtraExonData *eed2 = *((ExtraExonData **)b);

  return eed1->coords[0] - eed2->coords[0];
}

void ExtraExonData_free(ExtraExonData *eed) {
  free(eed->coords);
  free(eed);
}

typedef struct IntronCoordsStruct {
  long prevExonEnd;
  long nextExonStart;
  int strand;
  int isCanonical;
  int score;
} IntronCoords;

IntronCoords *IntronCoords_new(long prevExonEnd, long nextExonStart, int strand, int isCanon, int score) {
  IntronCoords *ic;

  if ((ic = (IntronCoords *)calloc(1,sizeof(IntronCoords))) == NULL) {
    fprintf(stderr,"Failed allocating IntronCoords\n");
    exit(EXIT_MEMORY);
  }
  ic->prevExonEnd   = prevExonEnd;
  ic->nextExonStart = nextExonStart;
  ic->strand        = strand;
  ic->isCanonical   = isCanon;
  ic->score         = score;

  return ic;
}

int IntronCoords_intronStartCompFunc(const void *a, const void *b) {
  IntronCoords *ic1 = *((IntronCoords **)a);
  IntronCoords *ic2 = *((IntronCoords **)b);

  return ic1->prevExonEnd - ic2->prevExonEnd;
}

void IntronCoords_free(IntronCoords *ic) {
  free(ic);
}

ModelCluster *ModelCluster_new() {
  ModelCluster *mc;

  if ((mc = (ModelCluster *)calloc(1,sizeof(ModelCluster))) == NULL) {
    fprintf(stderr,"Failed allocating ModelCluster\n");
    exit(EXIT_MEMORY);
  }

  return mc;
}

void ModelCluster_free(ModelCluster *mc) {
//  fprintf(stderr,"Freeing model cluster with ");
  if (mc->models) {
//    fprintf(stderr, "%d models", Vector_getNumElement(mc->models));
    Vector_setFreeFunc(mc->models, Model_free);
    Vector_free(mc->models);
  } else {
//    fprintf(stderr, "no models");
  }
// Final  
  if (mc->finalModels) {
//    fprintf(stderr, " and %d finalModels\n", Vector_getNumElement(mc->finalModels));
    //Vector_setFreeFunc(mc->finalModels, Gene_free);
    int i;
    for (i=0; i<Vector_getNumElement(mc->finalModels); i++) {
      Gene *g = Vector_getElementAt(mc->finalModels, i);
      Transcript *t = Gene_getTranscriptAt(g, 0);

      if (!(g->flags & RSGGENE_KEEP)) {
        if (t->extraData) TranscriptExtraData_free(t->extraData);
        Gene_free(g);
      }
    }
    Vector_free(mc->finalModels);
  } else {
//    fprintf(stderr, " and no finalModels\n");
  }
  free(mc);
}

SetFuncData *SetFuncData_new(void *func, int type) {
  SetFuncData *sfd;

  if ((sfd = (SetFuncData *)calloc(1,sizeof(SetFuncData))) == NULL) {
    fprintf(stderr,"Failed allocating SetFuncData\n");
    exit(EXIT_MEMORY);
  }
  sfd->setFunc.setIntValue = func;
  sfd->type = type;

  return sfd;
}


void SetFuncData_free(SetFuncData *sfd) {
  free(sfd);
}



#define RSG_DRIVER
#ifdef RSG_DRIVER

void RefineSolexaGenes_usage(int exit_code) {
  printf("RefineSolexaGenes \n"
         "  -c --config_file RefineSolexaGenes configuration file to read from\n"
         "  -i --input_id    Input id (slice name) to run on eg. chromosome:Oar_v3.1:17\n"
         "  -l --logic_name  Logic name for analysis block to run from configuration file\n"
         "  -d --dry_run     If specified, don't write to output db\n"
         "  -u --ucsc_naming If specified, add chr the name of sequence\n"
         "  -t --threads     Number of threads to use when reading the BAM files, default is 1\n"
         "  -v --verbosity   Verbosity level (int)\n"
         "  -V --version     Version of RefineSolexaGenes\n"
         "  -h --help        This help...\n"
         "\n"
         "Notes:\n"
         "  If no genes were written in the database, exit code is %i\n",
         NO_GENES_WRITTEN
         );
  exit(exit_code);
}

void RefineSolexaGenes_version() {
  printf("Version: %s\n", RSGVERSION);
  exit(EXIT_SUCCESS);
}

int main(int argc, char *argv[]) {
  if (argc == 1) RefineSolexaGenes_usage(EX_USAGE);
  initEnsC(argc, argv);

  logfp = stderr;


  char *logicName;
  char *configFile;
  char *inputId;
  int   dryRun     = 0;
  int   threads  = 1;
  int   verbosity  = 1;
  int   ucsc_naming = 0;

  int argNum = 1;
  while (argNum < argc) {
    char *arg = argv[argNum];
    char *val;

// Ones without a val go here
   if (!strcmp(arg, "-V") || !strcmp(arg,"--version")) {
     RefineSolexaGenes_version();
   }
   else if (!strcmp(arg, "-h") || !strcmp(arg,"--help")) {
     RefineSolexaGenes_usage(EXIT_SUCCESS);
   }
   else if (!strcmp(arg, "-d") || !strcmp(arg,"--dry_run")) {
      dryRun = 1;
   } else if (!strcmp(arg, "-u") || !strcmp(arg,"--ucsc_naming")) {
      ucsc_naming = 1;
    } else {
// Ones with a val go in this block
      if (argNum == argc-1) {
        fprintf(stderr, "Error: Expected a value after last command line argument\n");
        RefineSolexaGenes_usage(EX_USAGE);
      }

      val = argv[++argNum];
  //    printf("%s %s\n",arg,val);

      if (!strcmp(arg, "-c") || !strcmp(arg,"--config_file")) {
        StrUtil_copyString(&configFile,val,0);
      } else if (!strcmp(arg, "-i") || !strcmp(arg,"--input_id")) {
        StrUtil_copyString(&inputId,val,0);
      } else if (!strcmp(arg, "-l") || !strcmp(arg,"--logic_name")) {
        StrUtil_copyString(&logicName,val,0);
      } else if (!strcmp(arg, "-t") || !strcmp(arg,"--threads")) {
        threads = atoi(val);
      } else if (!strcmp(arg, "-v") || !strcmp(arg,"--verbosity")) {
        verbosity = atoi(val);
      } else {
        fprintf(stderr,"Error in command line at %s\n\n",arg);
        RefineSolexaGenes_usage(EX_USAGE);
      }
    }

    argNum++;
  }
  int exit_code = EXIT_SUCCESS;

  if (verbosity > 0) {
    printf("Program for generating transcript models from RNASeq read alignments\n"
           "Steve M.J. Searle.  searle@sanger.ac.uk  Last update Aug 2013.\n"
           "Thibaut Hourlier                         Last update %s\n", __DATE__);
  }

  RefineSolexaGenes *rsg = RefineSolexaGenes_new(configFile, logicName);
  if (verbosity > 0) RefineSolexaGenes_dumpConfig(rsg);
  RefineSolexaGenes_setInputId(rsg, inputId);
  RefineSolexaGenes_setDryRun(rsg, dryRun);
  RefineSolexaGenes_setThreads(rsg, threads);
  RefineSolexaGenes_setVerbosity(rsg, verbosity);
  RefineSolexaGenes_setUcscNaming(rsg, ucsc_naming);

//  double consLims[]    = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 20.0, 50.0 };
//  double nonConsLims[] = { 1.0, 5.0, 10.0, 15.0, 20.0, 50.0 };
//  double consLims[]    = { 1.0 };
//  double nonConsLims[] = { 5.0 };
//  double consLims[]    = { 5.0 };
//  double nonConsLims[] = { 5.0 };

  Vector *consLims = RefineSolexaGenes_getConsLims(rsg);
  Vector *nonConsLims = RefineSolexaGenes_getNonConsLims(rsg);
  double restartNonConsLim = RefineSolexaGenes_getRestartNonConsLim(rsg);

  // Check its a valid non cons lim to restart from (one that exists in nonConsLims)
  if (restartNonConsLim > -0.1) {
    int found = 0;
    int i;
    for (i=0; i<Vector_getNumElement(nonConsLims) && !found; i++) {
      double nonConsLim = *(double *)Vector_getElementAt(nonConsLims,i);
      if (nonConsLim == restartNonConsLim) {
        found = 1;
      }
    }
    if (!found) {
      fprintf(stderr, "Error: Didn't find restartNonConsLim %lf in nonConsLims\n",restartNonConsLim);
      exit(EX_DATAERR);
    }
// Don't write introns on a restart - they should already have been written
    RefineSolexaGenes_setWriteIntrons(rsg, 0);
  }
  

  //Vector *outputSets = Vector_new();

  int i;
  for (i=0; i<Vector_getNumElement(consLims); i++) {
    double consLim = *(double *)Vector_getElementAt(consLims,i);
    RefineSolexaGenes_setConsLim(rsg, consLim);
    int j;
    for (j=0; j<Vector_getNumElement(nonConsLims); j++) {
      double nonConsLim = *(double *)Vector_getElementAt(nonConsLims, j);
      RefineSolexaGenes_setNonConsLim(rsg, nonConsLim);

      // If we're doing a restart then check if at start NonConsLim and skip if not
      if (restartNonConsLim > -0.1 && i==0 && restartNonConsLim != nonConsLim) {
        if (verbosity > 0) fprintf(stderr,"Skipping consLim %lf nonConsLim %lf on restart\n", consLim, nonConsLim);
        continue;
      }

      char typeName[1024];
      char *typePref = RefineSolexaGenes_getTypePrefix(rsg);

      //sprintf(typeName, "%s_%sc%d_nc%d", RefineSolexaGenes_getBestScoreType(rsg), typePref, (int)consLim, (int)nonConsLim);
      //RefineSolexaGenes_setBestScoreType(rsg, typeName);

      //sprintf(typeName, "%s_%sc%d_nc%d", RefineSolexaGenes_getSingleExonModelType(rsg), typePref, (int)consLim, (int)nonConsLim);
      //RefineSolexaGenes_setSingleExonModelType(rsg, typeName);

      //sprintf(typeName, "%s_%sc%d_nc%d", logicName, typePref, (int)consLim, (int)nonConsLim);
      //RefineSolexaGenes_setAnalysis(rsg, RefineSolexaGenes_createAnalysisObject(rsg, typeName));
      RefineSolexaGenes_setAnalysis(rsg, RefineSolexaGenes_createAnalysisObject(rsg, logicName));

      RefineSolexaGenes_fetchInput(rsg);
      RefineSolexaGenes_run(rsg);
    
      fprintf(stderr,"Ended up with %d models to write\n", (RefineSolexaGenes_getOutput(rsg) ? Vector_getNumElement(RefineSolexaGenes_getOutput(rsg)) : 0) );

      if (verbosity > 0) fprintf(stderr, "cons lim %f non cons lim %f\n", consLim, nonConsLim);
      if (RefineSolexaGenes_getOutput(rsg)) {
        if (verbosity > 0) {
          dumpGenes(RefineSolexaGenes_getOutput(rsg), 1);
        }
        //fprintf(stderr,"malloc stats before write\n");
        //#ifdef HAVE_LIBTCMALLOC
        //tc_malloc_stats();
        //#endif
        if ( ! RefineSolexaGenes_isDryRun(rsg)) {
          RefineSolexaGenes_writeOutput(rsg);
        } else {
          fprintf(stderr,"DRY RUN mode - NOT writing genes to output db\n");
        }
#ifdef HAVE_LIBTCMALLOC
        if (verbosity > 0) tc_malloc_stats();
#endif
        ProcUtil_timeInfo("end of loop iter");
        if (verbosity > 0) fprintf(stderr,"Number of exon clone calls = %d\n",nExonClone);

        // HACK HACK HACK
        // For now only write introns for first iteration by setting WriteIntrons to 0 after it.
        // This reduces database load. We could load them in later if needed (they are just the 
        // set which is loaded from the BAM file and filtered using the cons and noncons filter params),
        // so its not a lot of compute to generate them, just a lot of database load and potentially
        // very large dna_align_feature tables.
        RefineSolexaGenes_setWriteIntrons(rsg, 0);
  
        int k;
        for (k=0;k<Vector_getNumElement(rsg->output);k++) {
          Gene *g = Vector_getElementAt(rsg->output, k);
          Transcript *t = Gene_getTranscriptAt(g, 0);
          if (t->extraData) {
            TranscriptExtraData_free(t->extraData);
          }
/*
          int m;
          for (m=0; m<Transcript_getExonCount(t); m++) {
            Exon *e = Transcript_getExonAt(t, m);
            Exon_clearSeqCacheString(e);
          }
*/
        }
        Vector_setFreeFunc(rsg->output, Gene_free);
        Vector_free(rsg->output);
      }
      else {
        exit_code = NO_GENES_WRITTEN;
      }
      //Vector_addElement(outputSets, rsg->output);
      rsg->output = NULL;
    }
  }

/*
  for (i=0;i<Vector_getNumElement(outputSets); i++) {
    Vector *set = Vector_getElementAt(outputSets, i);
    if (set != NULL) {
      //dumpGenes(set, 1);
    } else {
      fprintf(stderr,"Empty output set\n");
    }
  }
  //EcoString_getInfo(ecoSTable);
*/
#ifdef HAVE_LIBTCMALLOC
  if (verbosity > 0) tc_malloc_stats();
#endif
  return exit_code;
}
#endif

void RefineSolexaGenes_dumpOutput(RefineSolexaGenes *rsg) {
  Vector *output = RefineSolexaGenes_getOutput(rsg);

  Vector_sort(output, SeqFeature_startCompFunc);

  int i;
  for (i=0;i<Vector_getNumElement(output);i++) {
    Gene *gene = Vector_getElementAt(output, i);
    fprintf(stderr, "Gene\t%ld\t%ld\t%d\n", Gene_getStart(gene), Gene_getEnd(gene), Gene_getStrand(gene));
    int j;
    for (j=0; j<Gene_getTranscriptCount(gene); j++) {
      Transcript *trans = Gene_getTranscriptAt(gene, j);
      fprintf(stderr, "Transcript\t%ld\t%ld\t%d\n", Transcript_getStart(trans), Transcript_getEnd(trans), Transcript_getStrand(trans));
      int k;
      for (k=0;k<Transcript_getExonCount(trans);k++) {
        Exon *exon = Transcript_getExonAt(trans, k);
        fprintf(stderr, "Exon\t%ld\t%ld\t%d\n", Exon_getStart(exon), Exon_getEnd(exon), Exon_getStrand(exon));
      }
    }
    fprintf(stderr,"\n");
  }
}

void RefineSolexaGenes_initSetFuncs(RefineSolexaGenes *rsg) {
  rsg->funcHash = StringHash_new(STRINGHASH_SMALL);

  StringHash_add(rsg->funcHash, "OUTPUT_DB", SetFuncData_new(RefineSolexaGenes_setOutputDb, CONFIG_TYPE_STRING));
  StringHash_add(rsg->funcHash, "INTRON_DB", SetFuncData_new(RefineSolexaGenes_setIntronDb, CONFIG_TYPE_STRING));
  StringHash_add(rsg->funcHash, "MODEL_DB", SetFuncData_new(RefineSolexaGenes_setModelDb, CONFIG_TYPE_STRING));
  StringHash_add(rsg->funcHash, "WRITE_INTRONS", SetFuncData_new(RefineSolexaGenes_setWriteIntrons, CONFIG_TYPE_INT));
  StringHash_add(rsg->funcHash, "MAX_RECURSIONS", SetFuncData_new(RefineSolexaGenes_setMaxRecursions, CONFIG_TYPE_INT));

  SetFuncData *logicNameFuncData = SetFuncData_new(RefineSolexaGenes_setLogicNames, CONFIG_TYPE_ARRAY);
  logicNameFuncData->subType = CONFIG_TYPE_STRING;
  StringHash_add(rsg->funcHash, "LOGICNAME", logicNameFuncData);
  
  SetFuncData *intronBamFilesFuncData = SetFuncData_new(RefineSolexaGenes_setIntronBamFiles, CONFIG_TYPE_LIST);
  intronBamFilesFuncData->subType = CONFIG_TYPE_GROUP;
  intronBamFilesFuncData->subFunc = RefineSolexaGenes_parseIntronBamFilesConfig;
  StringHash_add(rsg->funcHash, "INTRON_BAM_FILES", intronBamFilesFuncData);

  StringHash_add(rsg->funcHash, "MODEL_LN", SetFuncData_new(RefineSolexaGenes_setModelLogicName, CONFIG_TYPE_STRING));
  StringHash_add(rsg->funcHash, "ALIGNMENT_BAM", SetFuncData_new(RefineSolexaGenes_setAlignmentBamFile, CONFIG_TYPE_STRING));
  StringHash_add(rsg->funcHash, "RETAINED_INTRON_PENALTY", SetFuncData_new(RefineSolexaGenes_setRetainedIntronPenalty, CONFIG_TYPE_FLOAT));
  StringHash_add(rsg->funcHash, "MIN_INTRON_SIZE", SetFuncData_new(RefineSolexaGenes_setMinIntronSize, CONFIG_TYPE_INT));
  StringHash_add(rsg->funcHash, "MAX_INTRON_SIZE", SetFuncData_new(RefineSolexaGenes_setMaxIntronSize, CONFIG_TYPE_INT));
  StringHash_add(rsg->funcHash, "INTRON_OVERLAP_THRESHOLD", SetFuncData_new(RefineSolexaGenes_setFilterOnOverlapThreshold, CONFIG_TYPE_INT));
  StringHash_add(rsg->funcHash, "IS_ONE_THRESHOLD", SetFuncData_new(RefineSolexaGenes_setIsOneThreshold, CONFIG_TYPE_INT));
  StringHash_add(rsg->funcHash, "SINGLE_EXON_MODEL", SetFuncData_new(RefineSolexaGenes_setSingleExonModelType, CONFIG_TYPE_STRING));
  StringHash_add(rsg->funcHash, "MIN_SINGLE_EXON", SetFuncData_new(RefineSolexaGenes_setMinSingleExonLength, CONFIG_TYPE_INT));
  StringHash_add(rsg->funcHash, "SINGLE_EXON_CDS", SetFuncData_new(RefineSolexaGenes_setMinSingleExonCDSPercLength, CONFIG_TYPE_FLOAT));
  StringHash_add(rsg->funcHash, "STRICT_INTERNAL_SPLICE_SITES", SetFuncData_new(RefineSolexaGenes_setStrictInternalSpliceSites, CONFIG_TYPE_INT));
  StringHash_add(rsg->funcHash, "STRICT_INTERNAL_END_EXON_SPLICE_SITES", SetFuncData_new(RefineSolexaGenes_setStrictInternalEndSpliceSites, CONFIG_TYPE_INT));
  StringHash_add(rsg->funcHash, "BEST_SCORE", SetFuncData_new(RefineSolexaGenes_setBestScoreType, CONFIG_TYPE_STRING));
  StringHash_add(rsg->funcHash, "OTHER_ISOFORMS", SetFuncData_new(RefineSolexaGenes_setOtherIsoformsType, CONFIG_TYPE_STRING));
  StringHash_add(rsg->funcHash, "OTHER_NUM", SetFuncData_new(RefineSolexaGenes_setOtherNum, CONFIG_TYPE_INT));
  StringHash_add(rsg->funcHash, "MAX_NUM", SetFuncData_new(RefineSolexaGenes_setMaxNum, CONFIG_TYPE_INT));
  StringHash_add(rsg->funcHash, "BAD_MODELS", SetFuncData_new(RefineSolexaGenes_setBadModelsType, CONFIG_TYPE_STRING));
  StringHash_add(rsg->funcHash, "TRIM_UTR", SetFuncData_new(RefineSolexaGenes_setTrimUTR, CONFIG_TYPE_INT));
  StringHash_add(rsg->funcHash, "MAX_3PRIME_EXONS", SetFuncData_new(RefineSolexaGenes_setMax3PrimeExons, CONFIG_TYPE_INT));
  StringHash_add(rsg->funcHash, "MAX_3PRIME_LENGTH", SetFuncData_new(RefineSolexaGenes_setMax3PrimeLength, CONFIG_TYPE_INT));
  StringHash_add(rsg->funcHash, "MAX_5PRIME_EXONS", SetFuncData_new(RefineSolexaGenes_setMax5PrimeExons, CONFIG_TYPE_INT));
  StringHash_add(rsg->funcHash, "MAX_5PRIME_LENGTH", SetFuncData_new(RefineSolexaGenes_setMax5PrimeLength, CONFIG_TYPE_INT));
  StringHash_add(rsg->funcHash, "REJECT_INTRON_CUTOFF", SetFuncData_new(RefineSolexaGenes_setRejectIntronCutoff, CONFIG_TYPE_FLOAT));
  StringHash_add(rsg->funcHash, "TYPE_PREFIX", SetFuncData_new(RefineSolexaGenes_setTypePrefix, CONFIG_TYPE_STRING));

  SetFuncData *consLimsFuncData = SetFuncData_new(RefineSolexaGenes_setConsLims, CONFIG_TYPE_ARRAY);
  consLimsFuncData->subType = CONFIG_TYPE_FLOAT;
  StringHash_add(rsg->funcHash, "CONSLIMS", consLimsFuncData);

  SetFuncData *nonConsLimsFuncData = SetFuncData_new(RefineSolexaGenes_setNonConsLims, CONFIG_TYPE_ARRAY);
  nonConsLimsFuncData->subType = CONFIG_TYPE_FLOAT;
  StringHash_add(rsg->funcHash, "NONCONSLIMS", nonConsLimsFuncData);

  StringHash_add(rsg->funcHash, "RESTART_NONCONSLIM", SetFuncData_new(RefineSolexaGenes_setRestartNonConsLim, CONFIG_TYPE_FLOAT));
}

void RefineSolexaGenes_dumpConfig(RefineSolexaGenes *rsg) {
  fprintf(stderr, "OUTPUT_DB\t\t%s\n", RefineSolexaGenes_getOutputDb(rsg));
  fprintf(stderr, "INTRON_DB\t\t%s\n", RefineSolexaGenes_getIntronDb(rsg));
  fprintf(stderr, "MODEL_DB\t\t%s\n", RefineSolexaGenes_getModelDb(rsg));
  fprintf(stderr, "WRITE_INTRONS\t\t%d\n", RefineSolexaGenes_writeIntrons(rsg));
  fprintf(stderr, "MAX_RECURSIONS\t\t%d\n", RefineSolexaGenes_getMaxRecursions(rsg));

  fprintf(stderr,"LOGICNAME\t\t["); 
  Vector *logicNames = RefineSolexaGenes_getLogicNames(rsg);
  if (logicNames) {
    int i;
    for (i=0; i<Vector_getNumElement(logicNames); i++) {
      char *lName = Vector_getElementAt(logicNames, i);
      if (i) fprintf(stderr,", ");
      fprintf(stderr,"'%s'",lName);
    }
  }
  fprintf(stderr,"]\n"); 
  
  fprintf(stderr,"INTRON_BAM_FILES {\n"); 
  Vector *intronBamFiles = RefineSolexaGenes_getIntronBamFiles(rsg);
  if (intronBamFiles) {
    int i;
    for (i=0; i<Vector_getNumElement(intronBamFiles); i++) {
      IntronBamConfig *ibc = Vector_getElementAt(intronBamFiles, i);
      fprintf(stderr," (\n  FILE %s\n  MIXED_BAM %d\n  DEPTH %d\n",ibc->fileName, ibc->mixedBam, ibc->depth);
      fprintf(stderr,"  GROUPNAME [ ");
      if (ibc->groupNames) {
        int j;
        for (j=0; j<Vector_getNumElement(ibc->groupNames); j++) {
          char *gName = Vector_getElementAt(ibc->groupNames, j);
          if (j) fprintf(stderr,", ");
          fprintf(stderr,"'%s'", gName);
        }
      }
      fprintf(stderr," ]\n");
      fprintf(stderr," )\n");
    }
  }
  fprintf(stderr,"}\n"); 

  fprintf(stderr, "MODEL_LN\t\t%s\n", RefineSolexaGenes_getModelLogicName(rsg));
  fprintf(stderr, "RETAINED_INTRON_PENALTY\t%lf\n", RefineSolexaGenes_getRetainedIntronPenalty(rsg));
  fprintf(stderr, "MIN_INTRON_SIZE\t\t%d\n", RefineSolexaGenes_getMinIntronSize(rsg));
  fprintf(stderr, "MAX_INTRON_SIZE\t\t%d\n", RefineSolexaGenes_getMaxIntronSize(rsg));
  fprintf(stderr, "SINGLE_EXON_MODEL\t\t%s\n", RefineSolexaGenes_getSingleExonModelType(rsg));
  fprintf(stderr, "MIN_SINGLE_EXON\t\t%d\n", RefineSolexaGenes_getMinSingleExonLength(rsg));
  fprintf(stderr, "SINGLE_EXON_CDS\t\t%lf\n", RefineSolexaGenes_getMinSingleExonCDSPercLength(rsg));
  fprintf(stderr, "STRICT_INTERNAL_SPLICE_SITES\t%d\n", RefineSolexaGenes_strictInternalSpliceSites(rsg));
  fprintf(stderr, "STRICT_INTERNAL_END_EXON_SPLICE_SITES\t%d\n", RefineSolexaGenes_strictInternalEndSpliceSites(rsg));
  fprintf(stderr, "BEST_SCORE\t\t%s\n", RefineSolexaGenes_getBestScoreType(rsg));
  fprintf(stderr, "OTHER_ISOFORMS\t\t%s\n", RefineSolexaGenes_getOtherIsoformsType(rsg));
  fprintf(stderr, "OTHER_NUM\t\t%d\n", RefineSolexaGenes_getOtherNum(rsg));
  fprintf(stderr, "MAX_NUM\t\t%d\n", RefineSolexaGenes_getMaxNum(rsg));
  fprintf(stderr, "BAD_MODELS\t\t%s\n", RefineSolexaGenes_getBadModelsType(rsg));
  fprintf(stderr, "TRIM_UTR\t\t%d\n", RefineSolexaGenes_trimUTR(rsg));
  fprintf(stderr, "MAX_3PRIME_EXONS\t\t%d\n", RefineSolexaGenes_getMax3PrimeExons(rsg));
  fprintf(stderr, "MAX_3PRIME_LENGTH\t\t%d\n", RefineSolexaGenes_getMax3PrimeLength(rsg));
  fprintf(stderr, "MAX_5PRIME_EXONS\t\t%d\n", RefineSolexaGenes_getMax5PrimeExons(rsg));
  fprintf(stderr, "MAX_5PRIME_LENGTH\t\t%d\n", RefineSolexaGenes_getMax5PrimeLength(rsg));
  fprintf(stderr, "REJECT_INTRON_CUTOFF\t\t%lf\n", RefineSolexaGenes_getRejectIntronCutoff(rsg));
  fprintf(stderr, "TYPE_PREFIX\t\t%s\n", RefineSolexaGenes_getTypePrefix(rsg));

  fprintf(stderr,"CONSLIMS\t\t["); 
  Vector *consLims = RefineSolexaGenes_getConsLims(rsg);
  if (consLims) {
    int i;
    for (i=0; i<Vector_getNumElement(consLims); i++) {
      double consLim = *(double *)Vector_getElementAt(consLims, i);
      if (i) fprintf(stderr,", ");
      fprintf(stderr,"%lf",consLim);
    }
  }
  fprintf(stderr,"]\n"); 

  fprintf(stderr,"NONCONSLIMS\t\t["); 
  Vector *nonConsLims = RefineSolexaGenes_getNonConsLims(rsg);
  if (nonConsLims) {
    int i;
    for (i=0; i<Vector_getNumElement(nonConsLims); i++) {
      double nonConsLim = *(double *)Vector_getElementAt(nonConsLims, i);
      if (i) fprintf(stderr,", ");
      fprintf(stderr,"%lf",nonConsLim);
    }
  }
  fprintf(stderr,"]\n"); 
  fprintf(stderr, "RESTART_NONCONSLIM\t\t%lf\n", RefineSolexaGenes_getRestartNonConsLim(rsg));
}

RefineSolexaGenes *RefineSolexaGenes_new(char *configFile, char *logicName) {
  RefineSolexaGenes *rsg;
  
  if ((rsg = calloc(1,sizeof(RefineSolexaGenes))) == NULL) {
    fprintf(stderr,"Failed allocating RefineSolexaGenes\n");
    exit(EXIT_MEMORY);
  }

  RefineSolexaGenes_initSetFuncs(rsg);

  rsg->adaptorAliasHash = StringHash_new(STRINGHASH_SMALL);

  if (configFile) {
    RunnableDB_readAndCheckConfig(rsg, configFile, "Config.REFINESOLEXAGENES_CONFIG_BY_LOGIC", logicName);
  } else {
    fprintf(stderr, "WARNING: Running without reading config (config reading not implemented)\n");
    exit(EX_DATAERR);
  }

//  $self->read_and_check_config($REFINESOLEXAGENES_CONFIG_BY_LOGIC);

  // Hard limit to the number of possible paths to explore
// SMJS Think maybe raise the start to 100000
  RefineSolexaGenes_setRecursiveLimit(rsg, RefineSolexaGenes_getMaxRecursions(rsg));

  // initialise intron feature cash
  // Doesn't seem to be used
  //my %feature_hash;
  //$self->feature_cash(\%feature_hash);

  return rsg;
}

void RunnableDB_readAndCheckConfig(RefineSolexaGenes *rsg, char *configFile, char *blockName, char *logicName) {
  config_t cfg;
  config_setting_t *cfgBlock;

  config_init(&cfg);

  /* Read the file. If there is an error, report it and exit. */
  if (!config_read_file(&cfg, configFile)) {
    fprintf(stderr, "%s:%d - %s\n", config_error_file(&cfg),
            config_error_line(&cfg), config_error_text(&cfg));
    config_destroy(&cfg);
    exit(EX_CONFIG);
  }

  // For now specify location of Databases config in config file and read it here
  config_setting_t *databasesFileSetting = config_lookup(&cfg, "Config.DATABASES_FILE");
  if (databasesFileSetting == NULL) {
    fprintf(stderr,"Missing config setting DATABASES_FILE\n");
    exit(EX_CONFIG);
  }
  const char *databasesFile = config_setting_get_string(databasesFileSetting);
  RunnableDB_readDatabaseConfig(rsg, (char *)databasesFile);

  cfgBlock = config_lookup(&cfg, blockName);
  if (cfgBlock == NULL) {
    fprintf(stderr,"Missing config block %s\n", blockName);
    exit(EX_CONFIG);
  }

  
// HACK: For now use passed in logicName rather than doing through analysis
//  Utilities_parseConfig(rsg, cfgBlock, Analysis_getLogicName(RefineSolexaGenes_getAnalysis(rsg)), 0 /*ignoreThrow*/);
  Utilities_parseConfig(rsg, cfgBlock, logicName, 0 /*ignoreThrow*/);
}

/*
=head2 fetch_sequence

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB
  Arg [2]   : string, name
  Arg [3]   : Bio::EnsEMBL::DBAdaptor
  Arg [4]   : arrayref of logic_name if sequence is to be masked
  Arg [5]   : Boolean for softmasking if sequence is to be softmasked
  Function  : gets sequence from specifed database
  Returntype: Bio::EnsEMBL::Slice
  Exceptions: none
  Example   :

=cut
*/
// This should go in RunnableDB
Slice *RefineSolexaGenes_fetchSequence(RefineSolexaGenes *rsg, char *name, DBAdaptor *db, Vector *repeatMaskTypes, int softMask) {
  if (!db) {
    db = RefineSolexaGenes_getDb(rsg);
  }

  if (!name) {
    name = RefineSolexaGenes_getInputId(rsg);
  }

  
  SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(db); 
  Slice *slice = SliceAdaptor_fetchByName(sa, name);

  if (slice == NULL) {
    fprintf(stderr,"Failed to fetch slice %s\n",name);
    exit(EX_IOERR);
  }

  fprintf(stderr, "slice name = %s\n", Slice_getName(slice));

  if (repeatMaskTypes != NULL && Vector_getNumElement(repeatMaskTypes)) {
    fprintf(stderr,"RepeatMasked sequence fetching not implemented yet - bye\n");
    exit(EX_SOFTWARE);
/* NIY
    my $sequence = $slice->get_repeatmasked_seq($repeat_masking, $soft_masking);
    $slice = $sequence
*/
  }

  return slice;
}

/*
=head2 fetch_input

    Title        :   fetch_input
    Usage        :   $self->fetch_input
    Returns      :   nothing
    Args         :   none

=cut
*/
void RefineSolexaGenes_fetchInput(RefineSolexaGenes *rsg) {
  int verbosity = RefineSolexaGenes_getVerbosity(rsg);

  if (verbosity > 0) fprintf(stderr, "Fetch input\n");
 // Not used?? char *logic = Analysis_getLogicName(RefineSolexaGenes_getAnalysis(rsg));
  // fetch adaptors and store them for use later
  // explicitly attach the ref db

  RefineSolexaGenes_setDb(rsg, BaseGeneBuild_getDbAdaptor(rsg, "REFERENCE_DB", 0, 0));

  DBAdaptor *db = RefineSolexaGenes_getDb(rsg);

  if (RefineSolexaGenes_getIntronDb(rsg) && RefineSolexaGenes_getIntronDb(rsg)[0] != '\0') {
    RefineSolexaGenes_setIntronSliceAdaptor(rsg, DBAdaptor_getSliceAdaptor(BaseGeneBuild_getDbAdaptor(rsg, RefineSolexaGenes_getIntronDb(rsg), 0, 0)));
  }
// Unused  $self->repeat_feature_adaptor($self->db->get_RepeatFeatureAdaptor);

  RefineSolexaGenes_setGeneSliceAdaptor(rsg, DBAdaptor_getSliceAdaptor(BaseGeneBuild_getDbAdaptor(rsg, RefineSolexaGenes_getModelDb(rsg), 0, 0)));

  // want a slice and a full chromsome to keep everything on in the same coords
  Slice *slice = RefineSolexaGenes_fetchSequence(rsg, RefineSolexaGenes_getInputId(rsg), NULL, NULL, 0);

  Slice *chrSlice = NULL;
  if (RefineSolexaGenes_getChrSlice(rsg) == NULL) {
      if (verbosity > 0) fprintf(stderr, "Got %d genes\n", __LINE__);
    chrSlice = SliceAdaptor_fetchByRegion(DBAdaptor_getSliceAdaptor(db), 
                                                 "toplevel", 
                                                 Slice_getSeqRegionName(slice),
                                                 POS_UNDEF,
                                                 POS_UNDEF,
                                                 1,
                                                 NULL,
                                                 0);
  
      if (verbosity > 0) fprintf(stderr, "Got %d genes\n", __LINE__);
    RefineSolexaGenes_setChrSlice(rsg, chrSlice);
      if (verbosity > 0) fprintf(stderr, "Got %d genes\n", __LINE__);
  } else {
      if (verbosity > 0) fprintf(stderr, "Got %d genes\n", __LINE__);
    chrSlice = RefineSolexaGenes_getChrSlice(rsg);
      if (verbosity > 0) fprintf(stderr, "Got %d genes\n", __LINE__);
  }
      if (verbosity > 0) fprintf(stderr, "Got %d genes\n", __LINE__);

  // we want to fetch and store all the rough transcripts at the beginning - speeds things up
  // also we want to take out tiny introns - merge them into longer structures
  if (RefineSolexaGenes_getPrelimGenes(rsg) == NULL) {
      if (verbosity > 0) fprintf(stderr, "Got %d genes\n", __LINE__);
    Slice *geneSlice =  SliceAdaptor_fetchByRegion(RefineSolexaGenes_getGeneSliceAdaptor(rsg),
                                                   "toplevel",
                                                   Slice_getSeqRegionName(slice),
                                                   Slice_getStart(slice),
                                                   Slice_getEnd(slice),
                                                   1,
                                                   NULL,
                                                   0);
      if (verbosity > 0) fprintf(stderr, "Got %d  %s genes\n", __LINE__, Slice_getName(geneSlice));
    Vector *genes;
      if (verbosity > 0) fprintf(stderr, "Got %d genes\n", __LINE__);
    char *modelLogicName = RefineSolexaGenes_getModelLogicName(rsg);
    if (modelLogicName != NULL) {
      if (verbosity > 0) fprintf(stderr, "Got %d genes\n", __LINE__);
      genes = Slice_getAllGenes(geneSlice, modelLogicName, NULL, 1, NULL, NULL);
      if (verbosity > 0) fprintf(stderr,"Got %d genes with logic name %s\n", Vector_getNumElement(genes), modelLogicName);
    } else {
      if (verbosity > 0) fprintf(stderr, "Got %d genes\n", __LINE__);
      genes = Slice_getAllGenes(geneSlice, NULL, NULL, 1, NULL, NULL);
      if (verbosity > 0) fprintf(stderr, "Got %d genes\n",Vector_getNumElement(genes));
    }
      if (verbosity > 0) fprintf(stderr, "Got %d genes\n", __LINE__);
  
    Vector *prelimGenes = Vector_new();
    int i;
      if (verbosity > 0) fprintf(stderr, "Got %d genes\n", __LINE__);
    for (i=0; i<Vector_getNumElement(genes); i++) {
      Gene *gene = Vector_getElementAt(genes, i);
  
      // put them on the chromosome
      gene = Gene_transfer(gene, chrSlice);
  
      // reject genes that are from a different slice that overlap our slice at the start or end
      // say the models has to be > 10% on the slice
      long os = Slice_getStart(slice);
      if (Gene_getStart(gene) > Slice_getStart(slice)) {
        os = Gene_getStart(gene);
      }
      long oe = Slice_getEnd(slice);
      if (Gene_getEnd(gene) < Slice_getEnd(slice)) {
        oe = Gene_getEnd(gene);
      }
  
      long overlap = oe - os + 1;
      
      double gc = ((double)overlap / (double)Gene_getLength(gene) * 1000.0) / 10.0;
      double sc =  ((double)overlap / (double)Slice_getLength(slice) * 1000.0) /10.0;
      if ( gc <= 10 && sc <= 10) {
        if (verbosity > 0) {
          fprintf(logfp, "Gene (start %ld end %ld) has %lf %% overlap with the slice\n"
                  "Slice (%s) has %lf %% overlap with the gene\n\tREJECTING\n",
                  Gene_getStart(gene), Gene_getEnd(gene), gc, Slice_getName(slice), sc);
        }
        continue;
      } 
      Vector_addElement(prelimGenes, gene);
    }
    if (verbosity > 0) fprintf(stderr, "Got %d genes after filtering boundary overlaps\n", Vector_getNumElement(prelimGenes));

    // Fake up a single exon gene for the entire region as test for improving models where rough models are very bitty
/* Test - didn't seem to work
    Transcript *trans = Transcript_new();
    Transcript_setSlice(trans, chrSlice);
    Transcript_setStableId(trans,"fake_trans");

    int j;
    for (j=Slice_getSeqRegionStart(slice); j<Slice_getSeqRegionEnd(slice); j+=10000) {
      Exon *exon = Exon_new();
      //Exon_setStart(exon, Slice_getSeqRegionStart(slice));
      //Exon_setEnd(exon, Slice_getSeqRegionEnd(slice));
      Exon_setStart(exon, j);
      Exon_setEnd(exon, j+9999);
      Exon_setPhase(exon, -1);
      Exon_setEndPhase(exon, -1);
      Exon_setStrand(exon, -1);
      Exon_setSlice(exon, chrSlice);
      Exon_setStableId(exon, "fake_exon");
  
      DNAAlignFeature *support = DNAAlignFeature_new();
      DNAAlignFeature_setStart(support, Exon_getStart(exon));
      DNAAlignFeature_setEnd(support, Exon_getEnd(exon));
      DNAAlignFeature_setStrand(support, Exon_getStrand(exon));
      DNAAlignFeature_setHitStrand(support, 1);
      DNAAlignFeature_setScore(support, 100.0);
      char tmpStr[1024];
      sprintf(tmpStr,"%ldM", Exon_getLength(exon));
      DNAAlignFeature_setCigarString(support, tmpStr);
      DNAAlignFeature_setSlice(support, chrSlice);
      DNAAlignFeature_setHitSeqName(support, "fake_support");
  
      Exon_addSupportingFeature(exon, support);
  
      Transcript_addExon(trans, exon, 0);
    }

    Gene *regGene = Gene_new();
    Gene_addTranscript(regGene, trans);
    Gene_setSlice(regGene, chrSlice);
    Gene_setStableId(regGene,"fake");

    fprintf(stderr,"Fake gene start %ld end %ld strand %d, trans start %ld end %ld strand %d\n",
            Gene_getStart(regGene), Gene_getEnd(regGene), Gene_getStrand(regGene),
            Transcript_getStart(trans), Transcript_getEnd(trans), Transcript_getStrand(trans));

    Vector_addElement(prelimGenes, regGene);
*/
  
    // determine strandedness ( including splitting merged genes )
    RefineSolexaGenes_setPrelimGenes(rsg, prelimGenes);
  }

  Vector *intronBamFiles = RefineSolexaGenes_getIntronBamFiles(rsg);

  if (Vector_getNumElement(intronBamFiles) > 0) {
    int i;
    char region_name[1024];
    if (RefineSolexaGenes_getUcscNaming(rsg) == 0) {
      sprintf(region_name,"%s", Slice_getSeqRegionName(slice));
    } else {
      sprintf(region_name,"chr%s", Slice_getSeqRegionName(slice));
    }
    for (i=0; i<Vector_getNumElement(intronBamFiles); i++) {
      IntronBamConfig *intronBamConf = Vector_getElementAt(intronBamFiles, i);
      char *intronFile = intronBamConf->fileName;
      char region[2048];
      int ref;
      int begRange;
      int endRange;
    
      htsFile *sam = hts_open(intronFile, "rb");
      if (sam == NULL) {
        fprintf(stderr, "Bam file %s not found\n", intronFile);
        exit(EX_NOINPUT);
      }
      if (verbosity > 0) fprintf(stderr,"Opened bam file %s\n", intronFile);

      hts_set_threads(sam, RefineSolexaGenes_getThreads(rsg));
      if (verbosity > 0) fprintf(stderr,"Setting number of threads to %d\n", RefineSolexaGenes_getThreads(rsg));
      hts_idx_t *idx;
      idx = sam_index_load(sam, intronFile); // load BAM index
      if (idx == 0) {
        fprintf(stderr, "BAM index file is not available.\n");
        exit(EX_NOINPUT);
      }
      if (verbosity > 0) fprintf(stderr,"Opened bam index for %s\n", intronFile);

      long long count = 0;
      bam_hdr_t *header = bam_hdr_init();
      header = bam_hdr_read(sam->fp.bgzf);
      ref = bam_name2id(header, region_name);
      if (ref < 0) {
        fprintf(stderr, "Invalid region %s\n", region_name);
        exit(EX_DATAERR);
      }
      sprintf(region,"%s:%ld-%ld", region_name, Slice_getSeqRegionStart(slice),
                                 Slice_getSeqRegionEnd(slice));

      if (hts_parse_reg(region, &begRange, &endRange) == NULL) {
        fprintf(stderr, "Could not parse %s\n", region);
        exit(EX_DATAERR);
      }
      if (verbosity > 0) fprintf(stderr,"Parsed region for region %s\n", region);

      // need to seamlessly merge here with the dna2simplefeatures code
      RefineSolexaGenes_bamToIntronFeatures(rsg, intronBamConf, sam, header, idx, ref, begRange, endRange);

      hts_idx_destroy(idx);
      bam_hdr_destroy(header);
      hts_close(sam);
    }
  } else {
    // pre fetch all the intron features
    RefineSolexaGenes_dnaToIntronFeatures(rsg, Slice_getStart(slice), Slice_getEnd(slice));
  }
}

void  RefineSolexaGenes_run(RefineSolexaGenes *rsg) {
//  fprintf(stderr, "RUN\n");
  RefineSolexaGenes_refineGenes(rsg);
//  fprintf(stderr, "FILTER\n");
  RefineSolexaGenes_filterGenes(rsg);
//  fprintf(stderr, "DONE\n");
}

int isGeneDuplicated(RefineSolexaGenes *rsg, Vector *indexes, Vector *genes, Gene *gene) {
  int i;
  int no_duplicate = 0;
  Gene *last_gene;
  Transcript *t;
  Transcript *lt;
  Exon *e;
  Exon *le;
  int num_exon;
  int lnum_exon;
  int j;
  for (i = 0; i< Vector_getNumElement(indexes); i++) {
    last_gene = Vector_getElementAt(genes, *(int *)Vector_getElementAt(indexes, i));
//    fprintf(stderr, "DEBUG: gene start %ld %ld %s\tLast: %ld %ld %s %d\n", Gene_getStart(gene), Gene_getEnd(gene), Gene_getBiotype(gene), Gene_getStart(last_gene), Gene_getEnd(last_gene), Gene_getBiotype(last_gene), *(int *)Vector_getElementAt(indexes, i));
    if (Gene_getStrand(gene) == Gene_getStrand(last_gene)
      && Gene_getStart(gene) == Gene_getStart(last_gene)
      && Gene_getEnd(gene) == Gene_getEnd(last_gene)) {
      t = Gene_getTranscriptAt(gene, 0);
      lt = Gene_getTranscriptAt(last_gene, 0);
      num_exon = Transcript_getExonCount(t);
      lnum_exon = Transcript_getExonCount(lt);
      if (num_exon == lnum_exon) {
        j = num_exon-1;
        do {
          e = Transcript_getExonAt(t, j);
          le = Transcript_getExonAt(lt, j);
          if (Exon_getStart(e) != Exon_getStart(le)
            || Exon_getEnd(e) != Exon_getEnd(le)) {
//            fprintf(stderr, "NOT DUPLICATE CHECK 2\n");
            ++no_duplicate;
            break;
          }
          --j;
        } while (j > 0);
        if (i > no_duplicate) {
//          fprintf(stderr, "DUPLICATE CHECK 2\n");
          if (strcmp(Gene_getBiotype(gene), RefineSolexaGenes_getBestScoreType(rsg)) == 0) {
//            fprintf(stderr, "UPDATE BIOTYPE FOR OLD GENE\n");
            updateDuplicatedGeneBiotype(last_gene);
            return *(int *)Vector_getElementAt(indexes, i);
          }
          else {
//            fprintf(stderr, "UPDATE BIOTYPE FOR GENE\n");
            updateDuplicatedGeneBiotype(gene);
          }
          return -1;
        }
      }
    }
  }
  return -2;
}

void updateDuplicatedGeneBiotype(Gene *gene) {
  if (!(gene->flags & RSG_DUPLICATE)) {
  //  int biotype_length = strlen(Gene_getBiotype(gene));
  //  char new_biotype[biotype_length+4];
  //  sprintf(new_biotype, "%s_dup", Gene_getBiotype(gene));
    gene->flags |= RSG_DUPLICATE;
    gene->flags ^= RSGGENE_KEEP;
  //  Gene_setBiotype(gene, new_biotype);
  }
}

void RefineSolexaGenes_filterGenes(RefineSolexaGenes *rsg) {
  int verbosity = RefineSolexaGenes_getVerbosity(rsg);
  Vector *genes = RefineSolexaGenes_getOutput(rsg);
  int i = 0;
  int to_delete_index = 0;
  int last_gene_index = 0;
  if (genes != NULL) {
    Vector_sort(genes, SeqFeature_startEndCompFunc);
    Gene *last_gene = Vector_getElementAt(genes, i);
    Transcript *t;
    Transcript *lt;
    int num_exon = 0;
    int lnum_exon = 0;
    int j;
    int duplicate;
    Exon *e;
    Exon *le;
    int biotype_length;
    Gene *gene;
    Vector *not_duplicate = Vector_new();
    Vector_addElement(not_duplicate, &i);
    for (i=1; i<Vector_getNumElement(genes); i++) {
      gene = Vector_getElementAt(genes, i);
  //    fprintf(stderr, "DEBUG: gene start %ld %ld %i %s\tLast: %ld %ld %i %s %d\n", Gene_getStart(gene), Gene_getEnd(gene), Gene_getStrand(gene), Gene_getBiotype(gene), Gene_getStart(last_gene), Gene_getEnd(last_gene), Gene_getStrand(last_gene), Gene_getBiotype(last_gene), last_gene_index);
      duplicate = -2;
      if (Gene_getStrand(gene) == Gene_getStrand(last_gene)
        && Gene_getStart(gene) == Gene_getStart(last_gene)
        && Gene_getEnd(gene) == Gene_getEnd(last_gene)) {
        t = Gene_getTranscriptAt(gene, 0);
        lt = Gene_getTranscriptAt(last_gene, 0);
        num_exon = Transcript_getExonCount(t);
        lnum_exon = Transcript_getExonCount(lt);
        if (num_exon == lnum_exon) {
          j = num_exon-1;
          duplicate = -1;
          do {
            e = Transcript_getExonAt(t, j);
            le = Transcript_getExonAt(lt, j);
            if (Exon_getStart(e) != Exon_getStart(le)
              || Exon_getEnd(e) != Exon_getEnd(le)) {
  //            fprintf(stderr, "NOT DUPLICATE\n");
              if (Vector_getNumElement(not_duplicate) > 1) {
  //              fprintf(stderr, "Checking other genes in the cluster\n");
                duplicate = isGeneDuplicated(rsg, not_duplicate, genes, gene);
              }
              else {
                duplicate = -2;
              }
              break;
            }
            --j;
          } while (j > 0);
  //        fprintf(stderr, "WE HAVE %i\n", duplicate);
          if (duplicate == -1) {
  //          fprintf(stderr, "DUPLICATE\n");
            if (strcmp(Gene_getBiotype(gene), RefineSolexaGenes_getBestScoreType(rsg)) == 0) {
              updateDuplicatedGeneBiotype(last_gene);
  //            fprintf(stderr, "REMOVING %s %ld %ld score %f index %i\n", Gene_getBiotype(last_gene), Gene_getStart(last_gene), Gene_getEnd(last_gene), Transcript_getScore(lt), last_gene_index);
  //            fprintf(stderr, "NEW IS BEST score %f index %i\t", Transcript_getScore(t), i);
              last_gene_index = i;
              last_gene = Vector_getElementAt(genes, last_gene_index);
  //            fprintf(stderr, "PROCESSED\n");
            }
            else {
  //            fprintf(stderr, "REMOVING %s %ld %ld score %f index %i\n", Gene_getBiotype(gene), Gene_getStart(gene), Gene_getEnd(gene), Transcript_getScore(t), i);
              updateDuplicatedGeneBiotype(gene);
            }
          }
          else if (duplicate > -1) {
  //          fprintf(stderr, "NEW BEST IS score %f index %i\t", Transcript_getScore(t), duplicate);
            last_gene_index = duplicate;
            last_gene = Vector_getElementAt(genes, last_gene_index);
  //          fprintf(stderr, "PROCESSED\n");
          }
          else {
  //          fprintf(stderr, "ADDING TO NOT DUPLICATE %i\t", i);
            Vector_addElement(not_duplicate, &i);
  //          fprintf(stderr, "PROCESSED\n");
          }
        }
        else {
              if (Vector_getNumElement(not_duplicate)) {
  //              fprintf(stderr, "Checking other genes in the cluster\n");
                duplicate = isGeneDuplicated(rsg, not_duplicate, genes, gene);
                if (duplicate > -1) {
  //                fprintf(stderr, "NEW BEST IS score %f index %i\t", Transcript_getScore(t), duplicate);
                  last_gene_index = duplicate;
                  last_gene = Vector_getElementAt(genes, last_gene_index);
  //                fprintf(stderr, "PROCESSED\n");
                }
                else {
  //                fprintf(stderr, "ADDING TO NOT DUPLICATE %i\t", i);
                  Vector_addElement(not_duplicate, &i);
  //                fprintf(stderr, "PROCESSED\n");
                }
              }
              else {
  //              fprintf(stderr, "not_duplicate has %i elements\n", Vector_getNumElement(not_duplicate));
              }
        }
      }
      else {
        last_gene_index = i;
        last_gene = Vector_getElementAt(genes, last_gene_index);
        Vector_free(not_duplicate);
        not_duplicate = Vector_new();
        Vector_addElement(not_duplicate, &i);
  //      fprintf(stderr, "BEST IS %s %ld %ld score %f index %i\n", Gene_getBiotype(gene), Gene_getStart(gene), Gene_getEnd(gene), Transcript_getScore(Gene_getTranscriptAt(gene, 0)), i);
      }
    }
    int k;
    int l;
    int min_exons = 0;
    Vector *clusters = Vector_new();
    int num_clusters = 0;
    int num_genes = 0;
    int model_added = 0;
    int m = 0;
    ModelCluster *cluster;
    for (i = 0; i < Vector_getNumElement(genes); i++) {
      gene = Vector_getElementAt(genes, i);
      if (gene->flags & RSG_DUPLICATE) {
        continue;
      }
      num_clusters = Vector_getNumElement(clusters);
      model_added = 0;
  //    fprintf(stderr, "DEBUG Index %i: %s %ld %ld %i\n", i, Gene_getBiotype(gene), Gene_getStart(gene), Gene_getEnd(gene), Gene_getStrand(gene));
      for (j = num_clusters-1; j > -1; j--) {
        cluster = Vector_getElementAt(clusters, j);
  //      fprintf(stderr, "DEBUG Cluster %i: %ld %ld %i\n", j, cluster->start, cluster->end, cluster->strand);
        num_genes = Vector_getNumElement(cluster->models);
        if (cluster->strand == Gene_getStrand(gene)) {
          if (cluster->start <= Gene_getEnd(gene) && cluster->end >= Gene_getStart(gene)) {
            // Gene overlaps, let's test the exons
            t = Gene_getTranscriptAt(gene, 0);
            for (k = 0; k < num_genes; k++) {
              last_gene = Vector_getElementAt(cluster->models, k);
              lt = Gene_getTranscriptAt(last_gene, 0);
              for (l = 0; l < Transcript_getExonCount(t); l++) {
                e = Transcript_getExonAt(t, l);
                for (m = 0; m < Transcript_getExonCount(lt); m++) {
                  le = Transcript_getExonAt(lt, m);
                  if (Exon_getEnd(e) < Exon_getStart(le)) {
                    break;
                  }
                  else if ((Exon_getStart(e) <= Exon_getEnd(le)) && (Exon_getEnd(e) >= Exon_getStart(le))) {
                    Vector_addElement(cluster->models, gene);
                    if (cluster->start > Gene_getStart(gene))
                      cluster->start = Gene_getStart(gene);
                    if (cluster->end < Gene_getEnd(gene))
                      cluster->end = Gene_getEnd(gene);
                    model_added = 1;
  //                  fprintf(stderr, "DEBUG ADDING %ld %ld %i to %i %ld %ld %i\n", Gene_getStart(gene), Gene_getEnd(gene), Gene_getStrand(gene), j, cluster->start, cluster->end, cluster->strand);
                    break;
                  }
                }
                if (model_added)
                  break;
              }
              if (model_added)
                break;
            }
            if (model_added)
              break;
          }
        }
      }
      if (!model_added) {
        cluster = ModelCluster_new();
        cluster->start = Gene_getStart(gene);
        cluster->end = Gene_getEnd(gene);
        cluster->strand = Gene_getStrand(gene);
        cluster->models = Vector_new();
        Vector_addElement(cluster->models, gene);
        Vector_addElement(clusters, cluster);
  //      fprintf(stderr, "DEBUG NEW %i %ld %ld %i\n", Vector_getNumElement(clusters)+1, cluster->start, cluster->end, cluster->strand);
      }
    }
  //  for (i = 0; i < Vector_getNumElement(clusters); i++) {
  //    cluster = Vector_getElementAt(clusters, i);
  //    fprintf(stderr, " Done\n");
  //    fprintf(stderr, "Cluster %i: %ld %ld %i\n", i, cluster->start, cluster->end, cluster->strand);
  //    for (j = 0; j < Vector_getNumElement(cluster->models); j++) {
  //      gene = Vector_getElementAt(cluster->models, j);
  //      fprintf(stderr, "Index %i: %ld %ld %i\n", j, Gene_getStart(gene), Gene_getEnd(gene), Gene_getStrand(gene));
  //    }
  //  }
    for (i = 0; i < Vector_getNumElement(clusters); i++) {
      cluster = Vector_getElementAt(clusters, i);
      Vector_sort(cluster->models, SeqFeat_lengthCompFunc);
  //    fprintf(stderr, "Num of models in cluster %i: %i\n", i, Vector_getNumElement(cluster->models));
      for (j = 0; j < Vector_getNumElement(cluster->models)-1; j++) {
        gene = Vector_getElementAt(cluster->models, j);
  //      fprintf(stderr, "WORKING ON %ld %ld %i\n", Gene_getStart(gene), Gene_getEnd(gene), Gene_getStrand(gene));
        for (k = j+1; k < Vector_getNumElement(cluster->models); k++) {
          t = Gene_getTranscriptAt(gene, 0);
          lt = Gene_getTranscriptAt((Gene *)Vector_getElementAt(cluster->models, k), 0);
          e = Transcript_getEndExon(t);
  //      fprintf(stderr, "%ld %ld ON %ld %ld %i\n", Exon_getStart(e), Exon_getEnd(e), Transcript_getStart(lt), Transcript_getEnd(lt), Transcript_getStrand(lt));
  //        fprintf(stderr, "last exon %i: %ld %ld %i\n", Transcript_getExonCount(t), Exon_getStart(e), Exon_getEnd(e), Exon_getStrand(e));
  //        fprintf(stderr, "Num of exon in model %i: %i\n", k, Transcript_getExonCount(lt));
          for (l = Transcript_getExonCount(lt)-2; l > -1; l--) {
            le = Transcript_getExonAt(lt, l);
            if (Exon_getStart(e) == Exon_getStart(le) && Exon_getEnd(e) >= Exon_getEnd(le)) {
              if(gene->flags & RSGGENE_KEEP) {
  //            fprintf(stderr, "END RETAINED %ld %ld %i\n", Gene_getStart(gene), Gene_getEnd(gene), Gene_getStrand(gene));
    //            biotype_length = strlen(Gene_getBiotype(gene));
    //            char new_biotype[biotype_length+4];
    //            sprintf(new_biotype, "%s_ret", Gene_getBiotype(gene));
    //            Gene_setBiotype(gene, new_biotype);
                gene->flags ^= RSGGENE_KEEP;
                // Set k to size of vector so for loop ends
                k = Vector_getNumElement(cluster->models);
                break;
              }
            }
            else if (Exon_getStart(e) > Exon_getStart(le)) {
  //            fprintf(stderr, "DONE %ld %ld %i\n", Exon_getStart(le), Exon_getEnd(le), Exon_getStrand(le));
              break;
            }
          }
          e = Transcript_getStartExon(t);
          for (l = 1; l < Transcript_getExonCount(lt)-1; l++) {
            le = Transcript_getExonAt(lt, l);
            if (Exon_getStart(e) <= Exon_getStart(le) && Exon_getEnd(e) == Exon_getEnd(le)) {
              if(gene->flags & RSGGENE_KEEP) {
  //              fprintf(stderr, "START RETAINED %ld %ld %i\n", Gene_getStart(gene), Gene_getEnd(gene), Gene_getStrand(gene));
    //            biotype_length = strlen(Gene_getBiotype(gene));
    //            char new_biotype[biotype_length+4];
    //            sprintf(new_biotype, "%s_ret", Gene_getBiotype(gene));
    //            Gene_setBiotype(gene, new_biotype);
                gene->flags ^= RSGGENE_KEEP;
                // Set k to size of vector so for loop ends
                k = Vector_getNumElement(cluster->models);
                break;
              }
            }
            else if (Exon_getEnd(e) < Exon_getEnd(le)) {
  //            fprintf(stderr, "DONE %ld %ld %i\n", Exon_getStart(le), Exon_getEnd(le), Exon_getStrand(le));
              break;
            }
          }
        }
      }
    }
    //A little bit hacky, I probably should have a clean function
    Vector *new_genes = Vector_new();
    for (i = 0; i < Vector_getNumElement(genes); i++) {
  //    RefineSolexaGenes_addToOutput(rsg, Vector_getElementAt(new_genes, i));
      gene = Vector_getElementAt(genes, i);
      if (gene->flags & RSGGENE_KEEP) {
        Vector_addElement(new_genes, Vector_getElementAt(genes, i));
  //      fprintf(stderr, "KEEPING %ld %ld %s\n", Gene_getStart((Gene *)Vector_getElementAt(genes, i)), Gene_getEnd((Gene *)Vector_getElementAt(genes, i)), Gene_getBiotype((Gene *)Vector_getElementAt(genes, i)));
      }
      else {
  //      fprintf(stderr, "DELETING %ld %ld %s\n", Gene_getStart((Gene *)Vector_getElementAt(genes, i)), Gene_getEnd((Gene *)Vector_getElementAt(genes, i)), Gene_getBiotype((Gene *)Vector_getElementAt(genes, i)));
        Gene_free(gene);
      }
    }
    Vector_free(genes);
    rsg->output = new_genes;
  }
}

/*
=head2 refine_genes

    Title        :   refine_genes
    Usage        :   $self->refine_genes
    Returns      :   nothing
    Args         :   none
    Description  :   Combines exons with introns in all possible combinations to 
                 :   Make a series of transcript models

=cut
*/

void RefineSolexaGenes_refineGenes(RefineSolexaGenes *rsg) {
  Vector *prelimGenes = RefineSolexaGenes_getPrelimGenes(rsg);
  int verbosity = RefineSolexaGenes_getVerbosity(rsg);

  int i;
// GENE:
  for (i=0; i<Vector_getNumElement(prelimGenes); i++) {

    if (verbosity > 1) fprintf(stderr,"STARTING NEW PRELIM GENE\n");
    Gene *gene = Vector_getElementAt(prelimGenes, i);
    Transcript *transcript = Gene_getTranscriptAt(gene, 0);
    int giveUpFlag = 0; // Note: This is used a long way down this routine - if set results in skipping to next gene which is why I put it here
   
    // hack taking out weeny models

    if (Transcript_getLength(transcript) < 300) {
      continue;
    }

    Vector *models = Vector_new();

    int singleExon = 0;
    // first run on the rev strand then on the fwd

//  STRAND: 
    int strand;
    for (strand = -1; strand<=1 && !giveUpFlag; strand+=2) {
// SMJS Maybe raise these defaults to 100000
      if ( RefineSolexaGenes_getRecursiveLimit(rsg) > RefineSolexaGenes_getMaxRecursions(rsg) ) {
        // set recursion to 10000 in case it was raised for a tricky gene
        RefineSolexaGenes_setRecursiveLimit(rsg, RefineSolexaGenes_getMaxRecursions(rsg));
        if (verbosity > 0) fprintf(stderr, "Warning: lowering recursive limit after complex gene\n");
      }

      if (verbosity > 1) fprintf(stderr,"Running on strand %d\n", strand);

// Doesn't seem to be used      StringHash *intronCount = StringHash_new(STRINGHASH_SMALL);
      Vector *exonIntron = Vector_new(); // A Vector of Vectors. Each Vector lists the possible introns for a particular exon
      
      StringHash *intronHash = StringHash_new(STRINGHASH_SMALL); // A hash containing all the introns, keyed on intron HitSeqName.
// Doesn't seem to be used      Vector *exonPrevIntron = Vector_new();
      StringHash *intronExon = StringHash_new(STRINGHASH_SMALL); // A Hash of Vectors. keyed on intron HitSeqName. Each element is a list of exons

      int mostRealIntrons = 0;
      double highestScore = 0;

      if (verbosity > 1) fprintf(stderr, "Gene %s : %ld %ld:\n", Gene_getStableId(gene), Gene_getStart(gene), Gene_getEnd(gene));

      Vector *exons = RefineSolexaGenes_mergeExons(rsg, gene, strand);
      Vector_sort(exons, SeqFeature_startCompFunc);

      int exonCount = Vector_getNumElement(exons);
// Doesn't seem to be used      Vector *fakeIntrons = Vector_new();
      StringHash *knownExons = StringHash_new(STRINGHASH_SMALL);

      long offset = 0;

      //fprintf(stderr,"exonCount = %d\n", exonCount);

//    EXON:   
      int j;
      for (j=0; j<exonCount; j++) {
        Exon *origExon = Vector_getElementAt(exons, j);
        Exon *exon = ExonUtils_cloneExon(origExon);
        int exonUsed = 0; // Has this cloned exon been added to exons (used for memory management)
 
        int retainedIntron = 0;
        int leftIntrons = 0;
        int rightIntrons = 0;

// Don't seem to be used so commented out
//        $exon->{'left_mask'} = 0;
//        $exon->{'right_mask'} = $exon->length;


        //fprintf(logfp, "Ex: %d : %ld %ld\n",j, Exon_getStart(exon), Exon_getEnd(exon));
        // make intron features by collapsing the dna_align_features

// Note in perl introns and offset are both returned - in C change to pass offset as pointer so can be modified
        Vector *introns = RefineSolexaGenes_fetchIntronFeatures(rsg, Exon_getSeqRegionStart(exon), Exon_getSeqRegionEnd(exon), &offset);

        Vector *leftConsIntrons = Vector_new();
        Vector *rightConsIntrons = Vector_new();
        Vector *leftNonConsIntrons = Vector_new();
        Vector *rightNonConsIntrons = Vector_new();
        Vector *filteredIntrons = Vector_new();
        Vector *retainedIntrons = Vector_new();

        int intronOverlap = 0;

//        fprintf(stderr, "Have %d introns\n", Vector_getNumElement(introns));
//      INTRON: 
        int k;
        for (k=0; k<Vector_getNumElement(introns); k++) {
          DNAAlignFeature *intron = Vector_getElementAt(introns, k);

          if (DNAAlignFeature_getStrand(intron) != strand) {
            //fprintf(stderr, "strand continue %d %d\n", DNAAlignFeature_getStrand(intron), strand);
            continue;
          } 
          if (DNAAlignFeature_getLength(intron) <= RefineSolexaGenes_getMinIntronSize(rsg)) {
            //fprintf(stderr, "min size continue\n");
            continue;
          } 
          if (DNAAlignFeature_getLength(intron) > RefineSolexaGenes_getMaxIntronSize(rsg)) {
            //fprintf(stderr, "max size continue\n");
            continue;
          } 

          // discard introns that splice over our exon
          if (DNAAlignFeature_getStart(intron) < Exon_getStart(exon) && 
              DNAAlignFeature_getEnd(intron) > Exon_getEnd(exon)) {
            intronOverlap++;
            continue;
          }
          //fprintf(stderr,"Got past continues\n");

          // check to see if this exon contains a retained intron
          if (DNAAlignFeature_getStart(intron) > Exon_getStart(exon) && 
              DNAAlignFeature_getEnd(intron) < Exon_getEnd(exon) ) {
          //fprintf(stderr,"in retained intron if \n");
            retainedIntron = 1;
// Hack hack hack

            Exon_addFlag(exon, RSGEXON_RETAINED);

            Vector_addElement(retainedIntrons, intron);
            //fprintf(stderr, "Added retained intron %s %d %d\n", DNAAlignFeature_getHitSeqName(intron), DNAAlignFeature_getStart(intron), DNAAlignFeature_getEnd(intron));
          } else {
            //fprintf(stderr, "Added normal intron %s %d %d\n", DNAAlignFeature_getHitSeqName(intron), DNAAlignFeature_getStart(intron), DNAAlignFeature_getEnd(intron));
          //fprintf(stderr,"in else \n");
            // we need to know how many consensus introns we have to the 
            // left and right in order to determine whether to put in 
            // a non consensus intron
            if (DNAAlignFeature_getEnd(intron) <= Exon_getEnd(exon)) {
              if (DNAAlignFeature_getFlags(intron) & RSGINTRON_NONCANON) {
                if (DNAAlignFeature_getScore(intron) > 1) {
                  Vector_addElement(leftNonConsIntrons, intron);
                }
              } else {
                Vector_addElement(leftConsIntrons, intron);
              }
            }
            if (DNAAlignFeature_getStart(intron) >= Exon_getStart(exon)) {
              if (DNAAlignFeature_getFlags(intron) & RSGINTRON_NONCANON) {
                if (DNAAlignFeature_getScore(intron) > 1) {
                  Vector_addElement(rightNonConsIntrons, intron);
                }
              } else {
                Vector_addElement(rightConsIntrons, intron);
              }
            }
          }
        }

/*
        fprintf(stderr, "Have %d left cons %d right cons %d left non cons %d right non cons\n", 
                Vector_getNumElement(leftConsIntrons),
                Vector_getNumElement(rightConsIntrons),
                Vector_getNumElement(leftNonConsIntrons),
                Vector_getNumElement(rightNonConsIntrons));
*/
        
        // Restrict internal exons splice sites to most common
        // that way our alt splices will all share the same boundaries
        // but have different combinations of exons
        if ( RefineSolexaGenes_strictInternalSpliceSites(rsg) || 
             // either we apply it equeally to all exons
             RefineSolexaGenes_strictInternalEndSpliceSites(rsg) ||
             // only apply to internal exons, leave out end exons
             ( !RefineSolexaGenes_strictInternalEndSpliceSites(rsg) &&
               ( Vector_getNumElement(leftConsIntrons)  + Vector_getNumElement(leftNonConsIntrons) ) > 0 &&
               ( Vector_getNumElement(rightConsIntrons) + Vector_getNumElement(rightNonConsIntrons)) > 0 )) {
          // pick best left splice
          long bestLeftSplice;
          double bestLeftScore = 0;

          int k;
          for (k=0; k<Vector_getNumElement(leftConsIntrons); k++) {
            DNAAlignFeature *intron = Vector_getElementAt(leftConsIntrons, k);

            if (bestLeftScore < DNAAlignFeature_getScore(intron)) {
              bestLeftScore = DNAAlignFeature_getScore(intron);
              bestLeftSplice = DNAAlignFeature_getEnd(intron);
            }
          }
          for (k=0; k<Vector_getNumElement(leftNonConsIntrons); k++) {
            DNAAlignFeature *intron = Vector_getElementAt(leftNonConsIntrons, k);

            if (bestLeftScore < DNAAlignFeature_getScore(intron)) {
              bestLeftScore = DNAAlignFeature_getScore(intron);
              bestLeftSplice = DNAAlignFeature_getEnd(intron);
            }
          }
          
          // pick best right  splice
          long bestRightSplice;
          double bestRightScore = 0;

          for (k=0; k<Vector_getNumElement(rightConsIntrons); k++) {
            DNAAlignFeature *intron = Vector_getElementAt(rightConsIntrons, k);

            if (bestRightScore < DNAAlignFeature_getScore(intron)) {
              bestRightScore = DNAAlignFeature_getScore(intron);
              bestRightSplice = DNAAlignFeature_getStart(intron);
            }
          }
          for (k=0; k<Vector_getNumElement(rightNonConsIntrons); k++) {
            DNAAlignFeature *intron = Vector_getElementAt(rightNonConsIntrons, k);

            if (bestRightScore < DNAAlignFeature_getScore(intron)) {
              bestRightScore = DNAAlignFeature_getScore(intron);
              bestRightSplice = DNAAlignFeature_getStart(intron);
            }
          }

          // filter out introns that pick other splice sites
          for (k=0; k<Vector_getNumElement(leftConsIntrons); k++) {
            DNAAlignFeature *intron = Vector_getElementAt(leftConsIntrons, k);

            if (DNAAlignFeature_getEnd(intron) == bestLeftSplice) {
              Vector_addElement(filteredIntrons, intron);
            }
          }
          for (k=0; k<Vector_getNumElement(leftNonConsIntrons); k++) {
            DNAAlignFeature *intron = Vector_getElementAt(leftNonConsIntrons, k);

            if (DNAAlignFeature_getEnd(intron) == bestLeftSplice) {
              Vector_addElement(filteredIntrons, intron);
            }
          }

          for (k=0; k<Vector_getNumElement(rightConsIntrons); k++) {
            DNAAlignFeature *intron = Vector_getElementAt(rightConsIntrons, k);

            if (DNAAlignFeature_getStart(intron) == bestRightSplice) {
              Vector_addElement(filteredIntrons, intron);
            }
          }
          for (k=0; k<Vector_getNumElement(rightNonConsIntrons); k++) {
            DNAAlignFeature *intron = Vector_getElementAt(rightNonConsIntrons, k);

            if (DNAAlignFeature_getStart(intron) == bestRightSplice) {
              Vector_addElement(filteredIntrons, intron);
            }
          }
        } else {
          
          // add non consensus introns only where there are no consensus introns
          Vector_append(filteredIntrons, leftConsIntrons);
          Vector_append(filteredIntrons, rightConsIntrons);

          if (Vector_getNumElement(leftConsIntrons) == 0) {
            Vector_append(filteredIntrons, leftNonConsIntrons);
          }
          if (Vector_getNumElement(rightConsIntrons) == 0) {
            Vector_append(filteredIntrons, rightNonConsIntrons);
          }
        }

        //fprintf(stderr, "Have %d exons %d filteredIntrons and %d retainedIntrons\n",  Vector_getNumElement(exons),
       //      Vector_getNumElement(filteredIntrons),
       //      Vector_getNumElement(retainedIntrons));
        
        // single exon models are a special case
        if ( Vector_getNumElement(exons) == 1 && 
             Vector_getNumElement(filteredIntrons) == 0 &&  
             Vector_getNumElement(retainedIntrons) == 0 ) {
          //# at least on this strand this model looks like a single exon
          singleExon += 1;
        }
        
        // we dont want to allow left and right introns to overlap - 
        // it leads to -ve length exons
        
        // we put all the retained introns in at the end we want to do all the 
        // entrances and exits to each exon before we look at whether its 
        // retained or not
        Vector_sort(retainedIntrons, SeqFeature_startCompFunc);

// Slight difference to perl - allocate a Vector for every exon in exonIntron. It will be empty if no introns for this exon, but easier than having a null pointer
        Vector *exIntj = NULL;

        // push @filtered_introns, @retained_introns;
//      INTRON:  
        for (k=0; k<Vector_getNumElement(filteredIntrons); k++) {
          DNAAlignFeature *intron = Vector_getElementAt(filteredIntrons, k);
          //# print STDERR "\t" . $intron->start . " " . $intron->end . " " . $intron->strand . " " . $intron->hseqname . " " . $intron->score . "\n";
          // becasue we make a new exons where we have a reatained intron to 
          // stop circular references we need to allow the final 
          // intron splicing out of the exon to be used more than once
          // by each new exon in fact

// Doesn't seem to be used          $intron_count{$intron->hseqname}++ unless $retained_intron;
//          StringHash_add(intronCount, DNAAlignFeature_getHitSeqName(intron),

          if (!StringHash_contains(intronHash, DNAAlignFeature_getHitSeqName(intron))) {
            StringHash_add(intronHash, DNAAlignFeature_getHitSeqName(intron), intron);
          }

          // only use each intron twice once at the end and once at the start of
          // an exon
          // exon_intron links exons to the intron on their right ignoring strand
          if (DNAAlignFeature_getEnd(intron) > Exon_getEnd(exon)) {
            if (exIntj == NULL) {
              exIntj = Vector_new();
              Vector_setElementAt(exonIntron, j, exIntj);
            }
            Vector_addElement(exIntj, intron);
          }

          // intron exon links introns to exons on their right ignoring strand
          if (DNAAlignFeature_getStart(intron) < Exon_getStart(exon)) {
            if (!StringHash_contains(intronExon, DNAAlignFeature_getHitSeqName(intron))) {
              Vector *ieVec = Vector_new();
              Vector_setFreeFunc(ieVec, free);
              StringHash_add(intronExon, DNAAlignFeature_getHitSeqName(intron), ieVec);
            }
            Vector *ieVec = StringHash_getValue(intronExon, DNAAlignFeature_getHitSeqName(intron));
            Vector_addElement(ieVec, long_new(j));
            // exon_prev_intron links exons to introns on their left ignoring strand
// Doesn't seem to be used push @{ $exon_prev_intron[$j]}  , $intron ;
          }
        }
        if (Vector_getNumElement(retainedIntrons) > 0) {
          //#print STDERR "Dealing with " . scalar( @retained_introns ) . " retained introns \n";
          //fprintf(stderr,  "Dealing with %d retained introns \n",Vector_getNumElement(retainedIntrons));
          Vector *newExons = Vector_new();
          exonUsed = 1;
          Vector_addElement(newExons, exon);

          // sort first by start then by end where start is the same
          Vector_sort(retainedIntrons, SeqFeature_startEndCompFunc);

          // now lets deal with any retained introns we have
//  RETAINED: 
          int m;
          for (m=0; m < Vector_getNumElement(retainedIntrons); m++) {
            DNAAlignFeature *intron = Vector_getElementAt(retainedIntrons, m);

            // we dont need to make all new exons for each alternate splice
            // check the intron is still retained given the new exons
            int retained = 1;
            int n;
            for (n=0; n<Vector_getNumElement(newExons); n++) {
              Exon *newExon = Vector_getElementAt(newExons, n);
              if  (DNAAlignFeature_getStart(intron) > Exon_getStart(newExon) && DNAAlignFeature_getEnd(intron) < Exon_getEnd(newExon)) {
              } else {
                retained = 0;
              }
              // Use an if instead of this - next RETAINED unless $retained;
            }

            if (retained) {
              double rejectScore = 0;
              // intron is within the exon - this is not a true exon but a retained intron
              if (DNAAlignFeature_getStart(intron) > Exon_getStart(exon) && 
                  DNAAlignFeature_getEnd(intron) < Exon_getEnd(exon) && 
                  DNAAlignFeature_getLength(intron) > RefineSolexaGenes_getMinIntronSize(rsg)) {
                // we are going to make a new exon and chop it up
                // add intron penalty
                //fprintf(stderr, "RETAINED INTRON PENALTY for %s before %f", DNAAlignFeature_getHitSeqName(intron), DNAAlignFeature_getScore(intron));
                rejectScore = DNAAlignFeature_getScore(intron) - RefineSolexaGenes_getRetainedIntronPenalty(rsg);
                // intron penalty is doubled for nc introns 
                //if (strstr(DNAAlignFeature_getHitSeqName(intron), "non canonical")) {
                if (DNAAlignFeature_getFlags(intron) & RSGINTRON_NONCANON) {
                  rejectScore = rejectScore - RefineSolexaGenes_getRetainedIntronPenalty(rsg);
                }
                //fprintf(stderr, " after %f\n",rejectScore);
                if (rejectScore < 1 ) {
                  // treat as single exon
                  if (Vector_getNumElement(exons) == 1 ) {
                    // at least on this strand this model looks like a single exon
                    singleExon += 1;
                  }
                } else {
                  //fprintf(stderr, "Exon %ld\t%ld has retained intron: %ld\t%ld\n",Exon_getStart(exon), Exon_getEnd(exon), DNAAlignFeature_getStart(intron), DNAAlignFeature_getEnd(intron));
                  // dont have circular references to exons or the paths
                  // will be infinite so clone this exon instead
                  // I guess we also want to keep the original exon too?
                  Exon *newExon1 = ExonUtils_cloneExon( exon );
                  Exon *newExon2 = ExonUtils_cloneExon( exon );
    
                  // chop it up a bit so it no longer overlaps the other introns
                  //fprintf(stderr, "TRIMMING EXON \n");
                  int length = DNAAlignFeature_getEnd(intron) - DNAAlignFeature_getStart(intron);
                  Exon_setEnd  (newExon1, DNAAlignFeature_getStart(intron) + ( length / 2 ) - 2);
                  Exon_setStart(newExon2, DNAAlignFeature_getEnd(intron) - ( length / 2 ) + 2);
                  
                  char keKey[1024];
                  sprintf(keKey,"%ld-%ld-%d",Exon_getStart(newExon1), Exon_getEnd(newExon1), Exon_getStrand(newExon1));
                  if (!StringHash_contains(knownExons, keKey)) {
                    Vector_addElement(newExons, newExon1);
                    StringHash_add(knownExons, keKey, &trueVal);
                  } else {
                    Exon_freeImpl(newExon1);
                  }
  
                  sprintf(keKey,"%ld-%ld-%d",Exon_getStart(newExon2), Exon_getEnd(newExon2), Exon_getStrand(newExon2));
                  if (!StringHash_contains(knownExons, keKey)) {
                    Vector_addElement(newExons, newExon2);
                    StringHash_add(knownExons, keKey, &trueVal);
                  } else {
                    Exon_freeImpl(newExon2);
                  }
                }
              }
            }
          }
          if (Vector_getNumElement(newExons) > 1) {
            // we want to split the score equally across the new exons
            
            for (k=0; k<Vector_getNumElement(newExons); k++) {
              Exon *e = Vector_getElementAt(newExons, k);
              int m;
              Vector *support = Exon_getAllSupportingFeatures(e);
              for (m=0; m<Vector_getNumElement(support); m++) { 
                BaseAlignFeature *d = Vector_getElementAt(support, m);
                BaseAlignFeature_setScore(d, BaseAlignFeature_getScore(d) / Vector_getNumElement(newExons));
              }
            }
            
            //splice( @exons,$i,1,@new_exons);
            Exon *removedExon = Vector_removeElementAt(exons, j);
            //Exon_freeImpl(removedExon);
            //fprintf(stderr,"Removing exon %p %ld %ld at %d\n", removedExon, Exon_getStart(removedExon), Exon_getEnd(removedExon), j);
            for (k=0; k<Vector_getNumElement(newExons); k++) {
              Exon *ex = Vector_getElementAt(newExons, k);
              Vector_insertElementAt(exons, j+k, ex);
              //fprintf(stderr,"Adding exon %p %ld %ld at %d\n", ex, Exon_getStart(ex), Exon_getEnd(ex), j+k);
            }
            //fprintf(logfp, "ADDED %d new exons\n", Vector_getNumElement(newExons));
            exonCount += (Vector_getNumElement(newExons) -1); // was $#new_exons;
            // make sure they are all stil sorted
            // This needs to sort reverse end for equal start to maintain the location of the retained intron version of
            // the exon before the cut ones
            Vector_sort(exons, SeqFeature_startRevEndCompFunc);
          }
          Vector_free(newExons);
        }
        Vector_free(introns);
        Vector_free(leftConsIntrons);
        Vector_free(rightConsIntrons);
        Vector_free(leftNonConsIntrons);
        Vector_free(rightNonConsIntrons);
        Vector_free(filteredIntrons);
        Vector_free(retainedIntrons);
        if (!exonUsed) Exon_freeImpl(exon);
      }
      StringHash_free(knownExons, NULL);
      
      // replaced with if - next unless @exon_intron;

      if (Vector_getNumElement(exonIntron) > 0) {

        // Loop around the path generation, 
        // if there are too many paths to process return undef
        // then re-run the path processing but with increasing strictness
        // where strictness = elimianating alternate low scoring introns
        if (verbosity > 1) fprintf(stderr, "STRAND %d BEFORE processPaths NUM EXONS %d num in exonIntron = %d\n", strand, Vector_getNumElement(exons), Vector_getNumElement(exonIntron));
        StringHash *paths = NULL;
        int strict = 0;
        while (paths == NULL) {
          paths = RefineSolexaGenes_processPaths(rsg, exons, exonIntron, intronExon, strict, &giveUpFlag );
  
          if (giveUpFlag) {
            //next GENE if $paths && $paths eq 'Give up';
            break; // Will use giveUpFlag to stop anything happening after here for this gene (except any necessary freeing of temp data)
          }
          strict++;
        }
  
        if (!giveUpFlag) {
          if (verbosity > 0) fprintf(stderr, "STRAND %d BEFORE COLLAPSING NUM PATHS = %d NUM EXONS %d\n", strand, StringHash_getNumValues(paths), Vector_getNumElement(exons));
    
          char **pathKeys = StringHash_getKeysNoCopy(paths);
//          char **pathKeys = StringHash_getKeys(paths);

          int nPath = StringHash_getNumValues(paths);


          // Note reason for making this array to store the paths to delete is to save space - the path strings can take up a lot of memory
          // so making a copy of the keys is expensive in both space and time terms. It's difficult to use the uncopied strings because we're
          // removing elements from the hash, so the key strings are freed, but could be used subsequently leading to a seg fault. By storing 
          // the ones to removed from the hash, and only removing them after the loop has finsihed executing we can avoid that problem.
          int nPathsToRemove = 0;
          char **pathsToRemove;
          if ((pathsToRemove = (char **)calloc(StringHash_getNumValues(paths), sizeof(char *))) == NULL) {
            fprintf(stderr, "Failed allocating pathsToRemove\n");
            exit(EXIT_MEMORY);
          }

          // Perl sorted them , not sure why - maybe longest first???
          // I didn't see any difference in results after commenting this out, so do without it (hopefully!)
          //qsort(pathKeys, nPath, sizeof(char *), StrUtil_stringCompFunc);

          int k;
          // lets collapse redundant paths
          for (k=0; k<nPath; k++) {
            char *path = pathKeys[k];

            char **array;
            int nArray;

            //fprintf(stderr,"PATHS %s\n", path);
            StrUtil_tokenizeByDelim(&array, &nArray, path, ".");

            char *startP;
            char *endP;
            char *middleP;
            
            char *start = (char *)calloc(100000,sizeof(char));
            char *end = (char *)calloc(100000,sizeof(char));
            char *middle = (char *)calloc(100000,sizeof(char));

            startP  = start;
            endP    = end;
            middleP = middle;

            int m;
            for (m=0; m<nArray; m++)  {
              //$start .= $array[m] ."."  unless k < 2;
              int lenArrayM = strlen(array[m]);
              
              
              //if (m>=2) sprintf(start,"%s%s.", start, array[m]);
              if (m>=2) {
                memcpy(startP, array[m], lenArrayM);
                startP+=lenArrayM;
                *startP++ = '.';
              }

              //$middle .= $array[m] ."." unless m < 2  or m >= $#array-1 ;
              //if (m>=2 && m<nArray-2) sprintf(middle,"%s%s.",middle, array[m]);
              if (m>=2 && m<nArray-2) {
                memcpy(middleP, array[m], lenArrayM);
                middleP+=lenArrayM;
                *middleP++ = '.';
              }

              //$end .= $array[m] . "." unless m >= $#array-1;
              //if (m<nArray-2) sprintf(end,"%s%s.", end, array[m]);
              if (m<nArray-2) {
                memcpy(endP, array[m], lenArrayM); // Note +1 is so we copy '\0' 
                endP+=lenArrayM;
                *endP++ = '.';
              }
            }
            *startP  = '\0';
            *endP    = '\0';
            *middleP = '\0';

            // remove redundancy from the array
            //delete $paths->{$start} if $start && $paths->{$start};
            char *keyP;
            if (start[0] != '\0' && (keyP = StringHash_getKey(paths, start)) != NULL) {
              //StringHash_remove(paths, start, NULL);
              //fprintf(stderr,"Adding start path %s\n", keyP);
              pathsToRemove[nPathsToRemove++] = keyP;
            }

            // delete $paths->{$end} if $end && $paths->{$end};
            if (end[0] != '\0' && (keyP = StringHash_getKey(paths, end)) != NULL) {
              //StringHash_remove(paths, end, NULL);
              //fprintf(stderr,"Adding end path %s\n", keyP);
              pathsToRemove[nPathsToRemove++] = keyP;
            }

            //delete $paths->{$middle} if $middle && $paths->{$middle};
            if (middle[0] != '\0' && (keyP = StringHash_getKey(paths, middle)) != NULL) {
              //StringHash_remove(paths, middle, NULL);
              //fprintf(stderr,"Adding middle path %s\n", keyP);
              pathsToRemove[nPathsToRemove++] = keyP;
            }

            for (m=0;m<nArray;m++) {
              free(array[m]);
            }
            free(array);
            free(start);
            free(end);
            free(middle);
          }
          
          //fprintf(stderr, "Have %d paths to remove\n", nPathsToRemove);
          qsort(pathsToRemove, nPathsToRemove, sizeof(char *), pointer_compFunc); // Note pointer sort NOT string sort
          for (k=0; k<nPathsToRemove; k++) {
            if (!k || pathsToRemove[k] != pathsToRemove[k-1]) {
              //fprintf(stderr, "Removing path |%s| (pointer %p)\n", pathsToRemove[k], pathsToRemove[k]);
              StringHash_remove(paths, pathsToRemove[k], NULL);
            }
          }

          free(pathKeys);
          free(pathsToRemove);
          
          if (verbosity > 0) fprintf(stderr, "AFTER COLLAPSING NUM PATHS  = %d NUM EXONS = %d\n", StringHash_getNumValues(paths), Vector_getNumElement(exons));
         
          Vector *strandModels = RefineSolexaGenes_makeModels(rsg, paths, strand, exons, gene, intronHash);
  
          Vector_append(models, strandModels);
          Vector_free(strandModels);
    
          //fprintf(stderr, "Now have %d models\n", Vector_getNumElement(models));
        }

        if (paths) {
          StringHash_free(paths, NULL);
          //fprintf(stderr,"Freeing paths %p\n",paths);
        }
      }

// NIY: Things to free
// exonIntron
      Vector_setFreeFunc(exonIntron, Vector_free);
      Vector_free(exonIntron);
// intronExon
      StringHash_free(intronExon, Vector_free);
// intronHash
      StringHash_free(intronHash, NULL);

      Vector_setFreeFunc(exons, Exon_freeImpl);
      Vector_free(exons);
    }

    //fprintf(stderr, "Give up flag = %d\n", giveUpFlag);

    // recursively recluster the models to identify 'other' models 
    // with no overlap to the 'best' model
    if (!giveUpFlag) {
      int modelCount = 0;
      Vector *clusteredModels;
      Vector *newClusters;
  
      clusteredModels = RefineSolexaGenes_reclusterModels(rsg, models, &newClusters);
  
      Vector *cleanClusters = Vector_new();
      
  // NIY: Nightmare to deal with memory freeing here
      if (newClusters != NULL && Vector_getNumElement(newClusters)) {
        while (Vector_getNumElement(newClusters)) {
          Vector_append(cleanClusters, clusteredModels);
          Vector_free(clusteredModels);
          Vector *oldNewClusters = newClusters;
          clusteredModels = RefineSolexaGenes_reclusterModels(rsg, newClusters, &newClusters);
          Vector_free(oldNewClusters);
          //fprintf(stderr, "Now have %d new clusters after reclustering\n", Vector_getNumElement(newClusters));
        }
        Vector_free(newClusters);
      } else if (newClusters != NULL) {
        Vector_free(newClusters);
      }
  
      if (clusteredModels != NULL) {
        if (Vector_getNumElement(clusteredModels)) {
          Vector_append(cleanClusters, clusteredModels);
        }
        Vector_free(clusteredModels);
      }
      
      // filter to identify 'best', 'other' and 'bad' models
      //fprintf(stderr,"XXXXXXXXXXXXX Have %d in cleanClusters\n", Vector_getNumElement(cleanClusters));
      //fprintf(stderr,"Number of final models in CLEAN clusters = %d\n", nFinal);

      RefineSolexaGenes_filterModels(rsg, cleanClusters);
      //Vector_setFreeFunc(cleanClusters, ModelCluster_free);
      Vector_free(cleanClusters);
  
      // process single exon models
      // if it has no introns on either strand
      if (RefineSolexaGenes_getSingleExonModelType(rsg) && singleExon == 2) {
        Vector *merged = RefineSolexaGenes_mergeExons(rsg, gene, 1);
        Exon *exon = Vector_getElementAt(merged, 0);

        Vector_free(merged);
  
        Transcript *singleExonModel = NULL;
        //fprintf(stderr, " Single exon = %d\n", singleExon);
    
        if (Exon_getLength(exon) + 40 >= RefineSolexaGenes_getMinSingleExonLength(rsg)) {
/*
          fprintf(stderr, "Passed length filter - exon start %ld end %ld length %ld\n", 
                  Exon_getStart(exon), Exon_getEnd(exon), Exon_getLength(exon));
*/
          // trim padding 
    // ?? Is padding always 20 ??
          Exon_setStart(exon,  Exon_getStart(exon) + 20);
          Exon_setEnd(exon,  Exon_getEnd(exon) - 20);
    
          // trim away strings of Ns from the start and  end
          // check start
          char *exSeq = Exon_getSeqString(exon);

          int nN;
          for (nN=0; nN<Exon_getLength(exon) && exSeq[nN] == 'N'; nN++) { }
          
          if (nN) {
            Exon_setStart(exon, Exon_getStart(exon) + nN);
          }

          // check end
          StrUtil_reverseString(exSeq, strlen(exSeq));
          for (nN=0; nN<Exon_getLength(exon) && exSeq[nN] == 'N'; nN++) { }
          if (nN) {
            Exon_setEnd(exon, Exon_getEnd(exon) - nN);
          }
          // NIY: Not sure if I need to free exSeq
         
          // get the cds
          Exon *fwdExon =  ExonUtils_cloneExon(exon);
          Exon_setStrand(fwdExon, 1);
    
          Exon *revExon = ExonUtils_cloneExon(exon);
    
          Transcript *fwdT = Transcript_new();
          Transcript_addExon(fwdT, fwdExon, 0);
          //my $fwd_t =  new Bio::EnsEMBL::Transcript(-EXONS => [$fwd_exon]);
          
// There were two more rounds of cloning here which seemed a bit excessive - I can't see why I need any 
       //   my $fwd_tran = TranscriptUtils_computeTranslation(clone_Transcript($fwdT));
          Transcript *fwdTran = TranslationUtils_computeTranslation(fwdT);
          
          long fwdTLen = 0;
          if (Transcript_getTranslation(fwdTran)) {
            fwdTLen = Translation_getGenomicEnd(Transcript_getTranslation(fwdTran)) - Translation_getGenomicStart(Transcript_getTranslation(fwdTran));
          }

          //fprintf(stderr, "FWD translation length %ld\n", fwdTLen);
          Transcript *revT = Transcript_new();
          Transcript_addExon(revT, revExon, 0);
          //my $rev_t =  new Bio::EnsEMBL::Transcript(-EXONS => [$rev_exon]);
          
         // my $rev_tran = compute_translation(clone_Transcript($rev_t));
          Transcript *revTran = TranslationUtils_computeTranslation(revT);
          
          long revTLen = 0;
          if (Transcript_getTranslation(revTran)) {
            revTLen = Translation_getGenomicEnd(Transcript_getTranslation(revTran)) - Translation_getGenomicStart(Transcript_getTranslation(revTran));
          }
    
          //fprintf(stderr, "REV translation length %ld\n", revTLen);
          if ( Transcript_getTranslation(fwdTran) &&  
               ( (double)fwdTLen / (double)Transcript_getLength(fwdTran)* 100.0 >= RefineSolexaGenes_getMinSingleExonCDSPercLength(rsg) &&
               fwdTLen >  revTLen )) {
            // keep this one
            singleExonModel =  fwdTran;

            Transcript_free(revTran);
            revTran = NULL;
          }
    
          if ( !singleExonModel && Transcript_getTranslation(revTran) &&  
               ( (double)revTLen / (double)Transcript_getLength(revTran)* 100.0 >= RefineSolexaGenes_getMinSingleExonCDSPercLength(rsg) &&
               revTLen >  fwdTLen )) {
            // keep this one
            singleExonModel =  revTran;

            Transcript_free(fwdTran);
            fwdTran = NULL;
          }

          if (singleExonModel == NULL) {
            if (fwdTran != NULL) Transcript_free(fwdTran);
            if (revTran != NULL) Transcript_free(revTran);
          }
    
          if (singleExonModel != NULL) {
            //fprintf(stderr, "Making single exon model\n");
            Transcript_setAnalysis(singleExonModel, RefineSolexaGenes_getAnalysis(rsg));
            Transcript_setVersion(singleExonModel, 1);
    
    // Was convert_to_genes which returned an array ref - make one which just returns a single gene
            Gene *newGene = TranscriptUtils_convertToGene(singleExonModel, Gene_getAnalysis(gene), NULL);
    
            Gene_setBiotype(newGene, RefineSolexaGenes_getSingleExonModelType(rsg));
    
            // score comes from exon supporting feature;
            Vector *support = Exon_getAllSupportingFeatures(exon);
            double score =  BaseAlignFeature_getScore((BaseAlignFeature *)Vector_getElementAt(support, 0));
            Exon_flushSupportingFeatures(exon);
    
            char stableId[2048];
            sprintf(stableId, "%s-v1-%d", Gene_getStableId(gene), (int)score);
            Gene_setStableId(newGene, stableId);
    
            //push @{$self->output} , $new_gene;
            RefineSolexaGenes_addToOutput(rsg, newGene);
          }
          Exon_free(exon);
        }
      }
    }
// NIY: Any freeing/tidying for a gene's worth of processing should go here

// NIY: Do we need to free models contents too??
    Vector_setFreeFunc(models, ModelCluster_free);
    //fprintf(stderr,"XXXXXXXXXXXXXXX Have %d in models\n", Vector_getNumElement(models));
    //fprintf(stderr,"Number of final models in all clusters = %d\n", nFinal);
    Vector_free(models);
#ifdef HAVE_LIBTCMALLOC
    MallocExtension_ReleaseFreeMemory();
#endif
  }
}

void RefineSolexaGenes_addToOutput(RefineSolexaGenes *rsg, Gene *gene) {
  if (rsg->output == NULL) {
    rsg->output = Vector_new();
  }
  gene->flags |= RSGGENE_KEEP;
  //fprintf(stderr,"Adding gene %s to output\n", Gene_getStableId(gene));
  Vector_addElement(rsg->output, gene);
}

Vector *RefineSolexaGenes_getOutput(RefineSolexaGenes *rsg) {
  return rsg->output;
}

// setAnalysis and getAnalysis should probably be in RunnableDB - put here for now
void RefineSolexaGenes_setAnalysis(RefineSolexaGenes *rsg, Analysis *analysis) {
  rsg->analysis = analysis;
}

Analysis *RefineSolexaGenes_getAnalysis(RefineSolexaGenes *rsg) {
  return rsg->analysis;
}

Analysis *RefineSolexaGenes_createAnalysisObject(RefineSolexaGenes *rsg, char *logicName) {
  DBAdaptor *outDb = BaseGeneBuild_getDbAdaptor(rsg, RefineSolexaGenes_getOutputDb(rsg), 0, 0);

  AnalysisAdaptor *aa = DBAdaptor_getAnalysisAdaptor(outDb);
  Analysis *analysis = AnalysisAdaptor_fetchByLogicName(aa, logicName);

  if (analysis != NULL) {
     //fprintf(stderr,"RETURNING existing analysis %p with logicName %s\n", analysis, logicName);
     return analysis;
  }
  // need to create analysis object first
  analysis = Analysis_new();
  char *chP;
  Analysis_setLogicName(analysis, StrUtil_copyString(&chP, logicName, 0)); 
  Analysis_setModule(analysis, StrUtil_copyString(&chP, "RefineSolexaGenes", 0));
  //fprintf(stderr,"MADE NEW analysis %p with logicName %s\n", analysis, logicName);

  return analysis;
}

/*
=head2 recluster_models

    Title        :   recluster_models
    Usage        :   $self->recluster_models(\@models);
    Returns      :   nothing
    Args         :   Array ref of array references of Bio::EnsEMBL::Transcript
    Description  :   reclusters 'other' models that have no overlap with 'best' models

=cut
*/
// called with 'models' which was immediately changed to clusters, so switched method arg to clusters
Vector *RefineSolexaGenes_reclusterModels(RefineSolexaGenes *rsg, Vector *clusters, Vector **retNewClusters) {
  Vector *newClusters   = Vector_new();
  Vector *finalClusters = Vector_new();

  //fprintf(stderr, "reclustering models - have %d clusters\n", Vector_getNumElement(clusters));

  int strand;
  for  (strand=1; strand >= -1; strand-=2) {
   // print "Running on strand $strand\n";
    
    int i;
    for (i=0; i<Vector_getNumElement(clusters); i++) {
      ModelCluster *cluster = Vector_getElementAt(clusters, i);

      if (cluster->finalModels == NULL) {
        continue;
      }

      Vector *strandedGenes = Vector_new();

      int j;
      for (j=0; j<Vector_getNumElement(cluster->finalModels); j++) {
        Gene *gene = Vector_getElementAt(cluster->finalModels, j);
        //fprintf(stderr, "GENE STRAND %d vs %d\n", Gene_getStrand(gene), strand);
        if (Gene_getStrand(gene) == strand) {
          Vector_addElement(strandedGenes, gene);
        }
      }

      if (Vector_getNumElement(strandedGenes) == 0) {
        Vector_free(strandedGenes);
        continue;
      }

      ModelCluster *strandedCluster = RefineSolexaGenes_recalculateCluster(rsg, strandedGenes);

/*
      fprintf(stderr, "Stranded Cluster %ld %ld strand %d with nfinal = %d\n",
              strandedCluster->start, strandedCluster->end, strandedCluster->strand, 
              Vector_getNumElement(strandedCluster->finalModels));
*/


      Gene *best = NULL;
      Vector *genes = Vector_new();
      Vector *otherGenes = Vector_new();

      long bestGenLen = 0;

      for (j=0; j<Vector_getNumElement(strandedCluster->finalModels); j++) {
        Gene *gene = Vector_getElementAt(strandedCluster->finalModels, j);

        if (Gene_getStrand(gene) == strand) {
          if (!strcmp(Gene_getBiotype(gene), RefineSolexaGenes_getBestScoreType(rsg))) {
            //$best =  $gene if $gene->biotype eq $self->BEST_SCORE;
// SMJS Note added some logic here to try to make consistent choice when multiple equal best (for perl/C comparisons)
//      Was previously just best = gene;
            if (best) {
              if (bestGenLen < Gene_getEnd(gene)-Gene_getStart(gene)+1) {
                best = gene;
                bestGenLen = Gene_getEnd(gene)-Gene_getStart(gene)+1;
              }
            } else {
              best = gene;
              bestGenLen = Gene_getEnd(gene)-Gene_getStart(gene)+1;
            }
            Vector_addElement(genes, gene);
          }
        }
      }

      if (best == NULL) {
        fprintf(stderr,"Error: No best model found\n");
        exit(EX_SOFTWARE);
      }

// Note: Moved these three lines out of loop below
      Transcript *bestTrans = Gene_getTranscriptAt(best, 0);
      Vector *bestExons = Vector_copy(Transcript_getAllExons(bestTrans));
// only need startComp but to make comparisons with perl need to use startEnd (to ensure same ordering between C and perl)
      Vector_sort(bestExons, SeqFeature_startEndCompFunc);

      // now recluster 
    //OTHERGENE: 
      for (j=0; j<Vector_getNumElement(strandedCluster->finalModels); j++) {
        Gene *gene = Vector_getElementAt(strandedCluster->finalModels, j);

/*
        fprintf(stderr,"Comparing gene %ld %ld %d (%s) to best gene %ld %ld %d (%s)\n",
                Gene_getStart(gene), Gene_getEnd(gene), Gene_getStrand(gene), Gene_getBiotype(gene),
                Gene_getStart(best), Gene_getEnd(best), Gene_getStrand(best), Gene_getBiotype(best)
               );
*/

        if (!strcmp(Gene_getBiotype(gene), RefineSolexaGenes_getBestScoreType(rsg))) {
          continue;
        }

        Transcript *geneTrans = Gene_getTranscriptAt(gene, 0);
        Vector *otherExons = Vector_copy(Transcript_getAllExons(geneTrans));
// only need startComp but to make comparisons with perl need to use startEnd (to ensure same ordering between C and perl)
        Vector_sort(otherExons, SeqFeature_startEndCompFunc);

        // exon overlap with a best model

        int doneBest = 0;
      //BESTEXON: 
        int k;
        for ( k=0; k < Vector_getNumElement(bestExons) && !doneBest; k++ ) {
          Exon *be = Vector_getElementAt(bestExons, k);

        //OTHEREXON: 
          int m;
          for (m=0 ; m < Vector_getNumElement(otherExons); m++) {
            Exon *oe = Vector_getElementAt(otherExons, m);

            if (Exon_getEnd(oe) < Exon_getStart(be)) {
              //next OTHEREXON;
              continue;
            }
// Big jump here - think about how to do
            if (Exon_getStart(oe) > Exon_getEnd(be)) {
              //next BESTEXON;
              break;
            } else {
              // does it have exon overlap with the best model
              if (Exon_getStart(be) <= Exon_getEnd(oe) && 
                  Exon_getEnd(be)   >= Exon_getStart(oe)) {
                // yes - store it and move on 
                Vector_addElement(genes, gene);
                
/*
                fprintf(stderr, "Overlap %ld %ld %d between exons %ld %ld %d and %ld %ld %d\n", 
                        Gene_getStart(gene), Gene_getEnd(gene), Gene_getStrand(gene),
                        Exon_getStart(be), Exon_getEnd(be), Exon_getStrand(be),
                        Exon_getStart(oe), Exon_getEnd(oe), Exon_getStrand(oe)
                       );
*/
  // Big jump here - think about how to do
                //next OTHERGENE;
                doneBest = 1;
                break;
              }
            }
          }
        }
        Vector_free(otherExons);
        // other model has no exon overlap with best model it needs to be in a new cluster
        if (!doneBest) {
          Vector_addElement(otherGenes, gene);
          //fprintf(stderr, "No overlap %ld %ld %d\n", Gene_getStart(gene), Gene_getEnd(gene), Gene_getStrand(gene));
        }
        //#print "No overlap " . $gene->start . " " , $gene->end . " " . $gene->strand ."\n";
      }
      Vector_free(bestExons);

      // now we need to fix the clusters
      if (Vector_getNumElement(otherGenes) > 0) {
        Vector_addElement(finalClusters, RefineSolexaGenes_recalculateCluster(rsg, genes));
        Vector_addElement(newClusters,   RefineSolexaGenes_recalculateCluster(rsg, otherGenes));
// NIY: Free strandedCluster???
// Try just the vectors
        //if (strandedCluster) {
           //Vector_free(strandedCluster->models);
           //Vector_free(strandedCluster->finalModels);
        //}
      } else {
        // keep it as it was
        Vector_free(genes);
        Vector_free(otherGenes);
        
        Vector_addElement(finalClusters, strandedCluster);
      }
// NIY: Free stuff??? - not sure this is the right place to do this!
    }
  }

  //fprintf(stderr,"CLUSTERS %d %d\n", Vector_getNumElement(finalClusters), Vector_getNumElement(newClusters));
  //return (\@final_clusters,\@new_clusters);
  *retNewClusters = newClusters;
  return finalClusters;
}

ModelCluster *RefineSolexaGenes_recalculateCluster(RefineSolexaGenes *rsg, Vector *genes) {
  ModelCluster *cluster = ModelCluster_new();
  Gene *firstGene = Vector_getElementAt(genes, 0);

  long start   = Gene_getStart(firstGene);
  long end     = Gene_getEnd(firstGene);
  long strand  = Gene_getStrand(firstGene);
  double score = 0;

  int i;
  for (i=0; i<Vector_getNumElement(genes); i++) {
    Gene *g = Vector_getElementAt(genes, i);
    if (Gene_getStart(g) < start) {
      start = Gene_getStart(g);
    }
    if (Gene_getEnd(g) > end) {
      end = Gene_getEnd(g);
    }

    Transcript *t = Gene_getTranscriptAt(g, 0);
    if (Transcript_getScore(t) > score) {
      score = Transcript_getScore(t);
    }
  }

  cluster->start  = start;
  cluster->end    = end;
  cluster->strand = strand;

  int nBest = 0;
  for (i=0; i<Vector_getNumElement(genes); i++) {
    Gene *g = Vector_getElementAt(genes, i);
    Transcript *t = Gene_getTranscriptAt(g, 0);

    if (Transcript_getScore(t) == score ) {
      Gene_setBiotype(g, RefineSolexaGenes_getBestScoreType(rsg));
      nBest++;
      //fprintf(stderr, "BEST SCORE %ld %ld %d %f\n", Gene_getStart(g), Gene_getEnd(g), Gene_getStrand(g), score);
    } else {
      Gene_setBiotype(g, RefineSolexaGenes_getOtherIsoformsType(rsg));
    }
  }

  cluster->finalModels = genes;

/*
  if (nBest > 1) {
    fprintf(stderr,"Ended up with multiple best scoring models in cluster (have %d)\n", nBest);
  }
*/

  return cluster;
}


/*
=head2 filter_models

    Title        :   filter_models
    Usage        :   $self->filter_models(\@models);
    Returns      :   nothing
    Args         :   Array ref of array references of Bio::EnsEMBL::Transcript
    Description  :   Labels or removes models overlapping better scoring models on the 
                     opposite strand

=cut
*/
// NOTE: Modifies output vector of RunnableDB
// Perl had models instead of clusters but then immediately assigned clusters = models, but then uses models for something else, so just call
// it clusters from start
void RefineSolexaGenes_filterModels(RefineSolexaGenes *rsg, Vector *clusters) {
  Vector *fwd = Vector_new();
  Vector *rev = Vector_new();

  int verbosity = RefineSolexaGenes_getVerbosity(rsg);

  Vector *models = Vector_new();

  //fprintf(stderr, "Filtering models - have %d clusters\n", Vector_getNumElement(clusters));
  int i;
  for (i=0; i<Vector_getNumElement(clusters); i++) {
    ModelCluster *cluster = Vector_getElementAt(clusters, i);

    if (cluster->finalModels == NULL) {
      if (verbosity > 0) fprintf(stderr,"Had a cluster with no final models\n");
      continue;
    }

    if (cluster->strand == 1) {
      Vector_addElement(fwd, cluster);
    }
    if (cluster->strand == -1) {
      Vector_addElement(rev, cluster);
    }
    Vector_addElement(models, cluster);
  }

  if (verbosity > 1) fprintf(stderr,"Have %d forward clusters and %d reverse clusters\n", Vector_getNumElement(fwd), Vector_getNumElement(rev));

  // overlaps
  for (i=0; i<Vector_getNumElement(fwd); i++) {
    ModelCluster *fc = Vector_getElementAt(fwd, i);
    int j;
    for (j=0; j<Vector_getNumElement(rev); j++) {
      ModelCluster *rc = Vector_getElementAt(rev, j);

      // one is within the other  or they are the same
      // they proably need to be rejected on the basis of coding overlap
      if (( fc->start >= rc->start && fc->end <= rc->end) || 
          ( rc->start >= fc->start && rc->end <= fc->end))  {
        
        // do they have coding overlap?
        Vector *fgs = fc->finalModels;
        Vector *rgs = rc->finalModels;

        //fprintf(stderr, "Overlapping clusters on two strands fwd one = %ld %ld (%d members) rev one = %ld %ld (%d members)\n", 
        //        fc->start, fc->end, Vector_getNumElement(fgs), rc->start, rc->end, Vector_getNumElement(rgs));

// Optimisation - precalculate translation lengths and translataebale exons for reverse genes so don't do millions of translate and translateable exon calls in loop!
        long *revLengths = NULL;
        if ((revLengths = (long *)calloc(Vector_getNumElement(rgs), sizeof(long))) == NULL) {
          fprintf(stderr, "Failed allocating array for reverse gene lengths\n");
          exit(EXIT_MEMORY);
        }
        Vector *revTranslateable = Vector_new();
        int k;
        for (k=0; k<Vector_getNumElement(rgs); k++) {
          Gene *rg = Vector_getElementAt(rgs, k);
          Transcript *rt = Gene_getTranscriptAt(rg, 0);

          if (!Transcript_getTranslation(rt)) {
            continue;
          }

/* Much slower than just summing exon lengths - only difference is getting seq removes stop codon so length can sometimes be 1 shorter.
          char *revTranslatedSeq = Transcript_translate(rt);
          revLengths[k] = strlen(revTranslatedSeq);
          free(revTranslatedSeq);
*/
          Vector *revTranslateableExons = Transcript_getAllTranslateableExons(rt);
          Vector_setElementAt(revTranslateable, k, revTranslateableExons);
// Do length by summing lengths of translateable exons - much quicker
          int m;
          for (m=0;m<Vector_getNumElement(revTranslateableExons);m++) {
            Exon *e = Vector_getElementAt(revTranslateableExons, m);
            revLengths[k] += Exon_getLength(e);
          }
          //fprintf(stderr,"rl[%d] = %d prior to /3\n", k, revLengths[k]);
          revLengths[k] /= 3;
// I added this after moving the code - not sure if its necessary
          if (revLengths[k] <= 100) {
            //fprintf(stderr,"bad RG len %ld %ld\n", Gene_getStart(rg), Gene_getEnd(rg));
            Gene_setBiotype(rg, "bad");
          } else {  
            //fprintf(stderr,"good RG len %ld for %ld %ld\n", revLengths[k], Gene_getStart(rg), Gene_getEnd(rg));
          }

        }

        // do they have coding overlap?
//      FG: 
        for (k=0; k<Vector_getNumElement(fgs); k++) {
          Gene *fg = Vector_getElementAt(fgs, k);
          Transcript *ft = Gene_getTranscriptAt(fg, 0);
          int done = 0;
          
          if (!Transcript_getTranslation(ft)) {
            if (verbosity > 0) fprintf(stderr,"!!!!!!!!!!!! No translation - continue\n");
            continue;
          }

/* Much slower than just summing exon lengths - only difference is getting seq removes stop codon so length can sometimes be 1 shorter.
          char *fwdTranslatedSeq = Transcript_translate(ft);
          int lenFwdTranslation = strlen(fwdTranslatedSeq);
          free(fwdTranslatedSeq);
*/

          Vector *ftTranslateableExons = Transcript_getAllTranslateableExons(ft);
          int lenFwdTranslation = 0;
          int m;
          for (m=0;m<Vector_getNumElement(ftTranslateableExons);m++) {
            Exon *e = Vector_getElementAt(ftTranslateableExons, m);
            lenFwdTranslation += Exon_getLength(e);
          }
          lenFwdTranslation /= 3;

          if (lenFwdTranslation <= 100) {
            Gene_setBiotype(fg, "bad");
            // Do an else instead next FG;
          } else {
//            RG: 
            int n;
            for (n=0; n<Vector_getNumElement(rgs) && !done; n++) {
              if (revLengths[n] == -1) {
                // Already set to bad - skip
                continue;
              }
              Gene *rg = Vector_getElementAt(rgs, n);
              Transcript *rt = Gene_getTranscriptAt(rg, 0);

// Next what?????
              if (!Transcript_getTranslation(rt)) {
                continue;
              }

              int lenRevTranslation = revLengths[n];

              if (lenRevTranslation <= 100) {
                Gene_setBiotype(rg, "bad");
                // Naughty - set revLengths to -1 to indicate we've flagged it as bad
                revLengths[n] = -1;
                  // Do an else instead next RG;
              } else {
                int m;
                for (m=0; m<Vector_getNumElement(ftTranslateableExons) && !done; m++) {
                  Exon *fe = Vector_getElementAt(ftTranslateableExons,  m);
                  int p;
                  Vector *rtTranslateableExons = Vector_getElementAt(revTranslateable, n);
                  for (p=0; p<Vector_getNumElement(rtTranslateableExons) && !done; p++) {
                    Exon *re = Vector_getElementAt(rtTranslateableExons, p);

                    if ( Exon_getStart(fe) <= Exon_getEnd(re) && 
                         Exon_getEnd(fe)  >=  Exon_getStart(re)) {
                      // coding overlap        

// Hack hack hack
                      if (Transcript_getScore(ft) < Transcript_getScore(rt)) {
                        //# get rid of / label the reverse genes 
                        //fprintf(stderr,"overlap - taking reverse - score rev %f v fwd %f\n", 
                        //        Transcript_getScore(rt), Transcript_getScore(ft));
                        Gene_setBiotype(fg, "bad");
                      } else {
                        //fprintf(stderr,"overlap - taking forward - score fwd %f v rev %f\n", 
                        //        Transcript_getScore(ft), Transcript_getScore(rt));
                        Gene_setBiotype(rg, "bad");
                      }
// Big jump - think about how to do
                      //next FG;
                      done = 1;
                    }
                  }
                  // Now Done after loop
                  //Transcript_freeAdjustedTranslateableExons(rt, rtTranslateableExons);
                  //Vector_free(rtTranslateableExons);
                }
              }
            }
          }
          Transcript_freeAdjustedTranslateableExons(ft, ftTranslateableExons);
          Vector_free(ftTranslateableExons);
        }
        for (k=0; k<Vector_getNumElement(rgs); k++) {
          Gene *rg = Vector_getElementAt(rgs, k);
          Transcript *rt = Gene_getTranscriptAt(rg, 0);

          if (!Transcript_getTranslation(rt)) {
            continue;
          }
          Vector *rtTranslateableExons = Vector_getElementAt(revTranslateable, k);
          if (rtTranslateableExons != NULL) {
            Transcript_freeAdjustedTranslateableExons(rt, rtTranslateableExons);
            Vector_free(rtTranslateableExons);
          } else {
        //    fprintf(stderr,"NULL rtTranslateableExons\n");
        //    exit(EX_SOFTWARE);
          }
        }
        Vector_free(revTranslateable);
        free(revLengths);
      }
    }
  }
  Vector_free(fwd);
  Vector_free(rev);

  fprintf(stderr, " got %d model clusters\n", Vector_getNumElement(models));

  for (i=0; i<Vector_getNumElement(models); i++) {
    ModelCluster *cluster = Vector_getElementAt(models, i);

    StringHash *exonUseHash = StringHash_new(STRINGHASH_SMALL);
    StringHash *exonPattern = StringHash_new(STRINGHASH_SMALL);
    IDHash *exonStarts      = IDHash_new(IDHASH_MEDIUM);
    IDHash *exonEnds        = IDHash_new(IDHASH_MEDIUM);

    int count = 0;
    long translationStart = 2000000000;
    long translationEnd   = 0;

    int j;
    for (j=0; j<Vector_getNumElement(cluster->finalModels); j++) {
      Gene *gene = Vector_getElementAt(cluster->finalModels, j);
      Transcript *transcript = Gene_getTranscriptAt(gene, 0);

     if (Transcript_getTranslation(transcript) != NULL) {
        if (Transcript_getCodingRegionStart(transcript) < translationStart) {
          translationStart = Transcript_getCodingRegionStart(transcript);
        }
        if (Transcript_getCodingRegionEnd(transcript) > translationEnd) {
          translationEnd = Transcript_getCodingRegionEnd(transcript);
        }
      }

      TranscriptExtraData *ted = Transcript_getExtraData(transcript);
      
      if (StringHash_contains(exonUseHash, ted->exonUse)) {
        Gene_setBiotype(gene, "bad");
      }

      StringHash_add(exonUseHash, ted->exonUse, &trueVal);

      //fprintf(stderr, "TRANSCRIPT %d Exon use %s biotype %s\n", ted->depth, ted->exonUse, Gene_getBiotype(gene));
      int es = 0;
      int ee = 0;
      char *pattern = (char *)calloc(65550,sizeof(char));
      pattern[0] = '\0';
      Transcript *trans = Gene_getTranscriptAt(gene, 0);
      int k;
      for (k=0; k<Transcript_getExonCount(trans); k++) {
        Exon *exon = Transcript_getExonAt(trans, k);

        sprintf(pattern,"%s%ld:%ld:",pattern,Exon_getStart(exon),Exon_getEnd(exon));

        if (IDHash_contains(exonStarts, Exon_getStart(exon))) {
          es++;
        } else {
          IDHash_add(exonStarts, Exon_getStart(exon), &trueVal);
        }
        if (IDHash_contains(exonEnds, Exon_getEnd(exon))) {
          ee++;
        } else {
          IDHash_add(exonEnds, Exon_getEnd(exon), &trueVal);
        }
      }
      if ( ee == Transcript_getExonCount(trans) &&
           es == Transcript_getExonCount(trans)) { 
        // seen it before - or something very much like it
        Gene_setBiotype(gene, "bad");
        //fprintf(stderr, "CALLING it bad\n");
      }

      if (StringHash_contains(exonPattern, pattern)) {
        // seen it before - or something very much like it
        Gene_setBiotype(gene, "duplicate");
        //fprintf(stderr, "CALLING it duplicate\n");
      } else {
        StringHash_add(exonPattern, pattern, &trueVal);
      }

      free(pattern);
    } 
    IDHash_free(exonStarts, NULL);
    IDHash_free(exonEnds, NULL);
    StringHash_free(exonPattern, NULL);
    StringHash_free(exonUseHash, NULL);

    // promote "bad" models that have a cds as long as the best cds to 
    // alt isoforms
    Vector *finalModels = cluster->finalModels;
    int bestCds = 0;
    int g;
    for (g=0; g<Vector_getNumElement(finalModels); g++) {
      Gene *gene = Vector_getElementAt(finalModels, g);


      Transcript *transcript = Gene_getTranscriptAt(gene, 0);
  
/*
      fprintf(logfp, "%d - %.4f tran length %d\n", 
              g, Transcript_getScore(transcript), Transcript_getcDNACodingEnd(transcript) - Transcript_getcDNACodingStart(transcript) + 1);
*/
/*
      if (Transcript_getTranslation(transcript)) {
        fprintf(logfp,"cds         seq %s\n", Transcript_getTranslateableSeq(transcript));
        fprintf(logfp,"translation seq %s\n", Transcript_translate(transcript));
      }
*/
        //int n;
        //for (n=0;n<Transcript_getExonCount(transcript);n++) {
         // Exon *e = Transcript_getExonAt(transcript, n);
          //fprintf(stderr," Exon %ld %ld %d\n", Exon_getStart(e), Exon_getEnd(e), Exon_getStrand(e));
        //}

      if (g == 0) {
        // best scoring model 
        if  (Transcript_getTranslation(transcript) != NULL) {
          bestCds = Transcript_getcDNACodingEnd(transcript) - Transcript_getcDNACodingStart(transcript);
        }
      }
      if (Transcript_getTranslation(transcript) != NULL) {
        if (!strcmp(Gene_getBiotype(gene), "bad") && 
             Transcript_getCodingRegionStart(transcript) == translationStart &&
             Transcript_getCodingRegionEnd(transcript) == translationEnd ) {
          Gene_setBiotype(gene, RefineSolexaGenes_getOtherIsoformsType(rsg));
        }
        if (!strcmp(Gene_getBiotype(gene), "bad") && 
            Transcript_getcDNACodingEnd(transcript) - Transcript_getcDNACodingStart(transcript) > bestCds ) {
          Gene_setBiotype(gene, RefineSolexaGenes_getOtherIsoformsType(rsg));
        }
      } 
      if (!strcmp(Gene_getBiotype(gene), "bad")) {
        // change type to  a bad model if it is bad 
        // store it on output array if the bad type is defined
        if (RefineSolexaGenes_getBadModelsType(rsg) != NULL) {
          Gene_setBiotype(gene, RefineSolexaGenes_getBadModelsType(rsg));
          if (count <= RefineSolexaGenes_getOtherNum(rsg)) {
            RefineSolexaGenes_addToOutput(rsg, gene);
            //Vector_setElementAt(finalModels, g, NULL);
          }
        } else { // Free it
          //Gene_free(gene);
        }
      } else {
        // if NOT duplicate 
        if (strcmp(Gene_getBiotype(gene), "duplicate")) {
          if (!strcmp(Gene_getBiotype(gene), RefineSolexaGenes_getBestScoreType(rsg))) {
            // trim the UTR
            if (RefineSolexaGenes_pruneUTR(rsg, gene) == 0) {
              RefineSolexaGenes_addToOutput(rsg, gene);
            }
//            else {
//              Gene_free(gene);
//            }
            //Vector_setElementAt(finalModels, g, NULL);
          } else {
// Note here checking for whether other isoforms type is null, other places assume its not null
            if (RefineSolexaGenes_getOtherNum(rsg) && RefineSolexaGenes_getOtherIsoformsType(rsg) != NULL && strlen(RefineSolexaGenes_getOtherIsoformsType(rsg)) != 0 && count <= RefineSolexaGenes_getOtherNum(rsg)) {
              // trim the UTR
              if (RefineSolexaGenes_pruneUTR(rsg, gene) == 0) {
                RefineSolexaGenes_addToOutput(rsg, gene);
              }
//            else {
//              Gene_free(gene);
//            }
              //Vector_setElementAt(finalModels, g, NULL);
            } else { // Free it
              //Gene_free(gene);
            }
          }
        } else { // Free it
          //Gene_free(gene);
        }
      }
      if (!strcmp(Gene_getBiotype(gene), RefineSolexaGenes_getOtherIsoformsType(rsg)) || 
          !strcmp(Gene_getBiotype(gene), RefineSolexaGenes_getBestScoreType(rsg))) {
        count++;
      }
    }
  }
  Vector_free(models);
}


/*
=head2 make_models

    Title        :   
    Usage        :   $self->make_models($paths, $strand ,$exons,$gene, $intron_hash);
    Returns      :   Array ref of array references of Bio::EnsEMBL::Transcript
    Description  :   Turns abstract paths into Bio::EnsEMBL::Gene models. Paths are
                     clustered and sorted by score - only the top X models for 
                     each cluster of paths get built ( X is defined in config )

=cut
*/

Vector *RefineSolexaGenes_makeModels(RefineSolexaGenes *rsg, StringHash *paths, int strand, Vector *exons, Gene *gene, StringHash *intronHash) {
  // paths are stored as text - turn them into arrays of features "models"
  Vector *clusters = Vector_new();
  Vector *models = Vector_new();
// Unused Vector *genes = Vector_new();

  char **pathKeys = StringHash_getKeys(paths);
  int nPath = StringHash_getNumValues(paths);
  int verbosity = RefineSolexaGenes_getVerbosity(rsg);

  int i;
  for (i=0; i<nPath; i++) {
    char *path = pathKeys[i];

    char exonUse[8000]; exonUse[0] = '\0';

    Model *model;
    if ((model = (Model *)calloc(1,sizeof(Model))) == NULL) {
      fprintf(stderr,"Failed allocating model\n");
      exit(EXIT_MEMORY);
    } 

    model->features = Vector_new();
    double exonScore = 0;
    double intronScore = 0;

    char **featureStrings;
    int nFeatureString;

    StrUtil_tokenizeByDelim(&featureStrings, &nFeatureString, path, ".");

    //fprintf(stderr, "path to make model from is %s\n", path);

    int j;
    for (j=0; j<nFeatureString; j++)  {
      char *featureString = featureStrings[j];
      if (strstr(featureString, "canon")) { // Intron
        DNAAlignFeature *intron = StringHash_getValue(intronHash, featureString);
        if (intron == NULL) {
          fprintf(stderr, "Couldn't find intron in intronHash. Key = %s\n",featureString);
        }
        Vector_addElement(model->features, intron);
        intronScore += DNAAlignFeature_getScore(intron);
      } else { // Exon
        sprintf(exonUse, "%s%s.",exonUse, featureString);
        int exonInd;
        sscanf(featureString, "%d", &exonInd);
        Exon *exon = Vector_getElementAt(exons, exonInd);
        Vector_addElement(model->features, exon);

        Vector *support = Exon_getAllSupportingFeatures(exon);
        int k;
        for (k=0; k<Vector_getNumElement(support); k++) {
          DNAAlignFeature *daf = Vector_getElementAt(support, k);
          exonScore += DNAAlignFeature_getScore(daf);
        }
      }
    }

    for (j=0; j<nFeatureString; j++)  {
      free(featureStrings[j]);
    }
    free(featureStrings);


    double totalScore = ((int)exonScore)/100 + intronScore;
    // last elements are the strand and score
    
    StrUtil_copyString(&model->exonUse, exonUse, 0);
    model->totalScore = totalScore;
    Vector_addElement(models, model);

    // NIY: Free stuff
  }

  for (i=0; i<nPath; i++) {
    free(pathKeys[i]);
  }
  free(pathKeys);

  //fprintf(logfp, "Starting model_cluster\n");
  // now lets cluster the models so that they are non overlapping
  // and return the clusters arranged by score

  Vector *modelClusters = RefineSolexaGenes_makeModelClusters(rsg, models, strand);
  Vector_free(models);
  
  //fprintf(logfp, "Starting gene cycle\n");
  // Now we cycle through all the models and turn them into genes
  // we start with the highest scoring models and work backwards
  // until we have enough 
  int clusterCount = 0;
  for (i=0; i<Vector_getNumElement(modelClusters); i++) {
    ModelCluster *cluster = Vector_getElementAt(modelClusters, i);
    clusterCount++;
    Vector *trans = Vector_new();
    int version = 0;

    int strand = cluster->strand;
    // we want the array in the reverse order, highest scoring first
    // the score is the last array element 
    Vector *modelsByScore = Vector_copy(cluster->models);
    Vector_sort(modelsByScore, Model_reverseScoreCompFunc);

    //my @models_by_score =sort {$b->[-1] <=> $a->[-1]}  @{$cluster->{'models'}};
    
    // all the models with a particular score highest first
//  MODEL:   
    int j;
    for (j=0; j<Vector_getNumElement(modelsByScore); j++) {
      Model *model = Vector_getElementAt(modelsByScore, j);

      Vector *ises = Vector_new();
      // the score is the last array element 
      double s = model->totalScore;

      // list of the rough exons used in the model
      char *exonUse = model->exonUse;

      Vector *introns = Vector_new();

      int intronCount = 0;
      int nonConIntrons = 0;
      double intronScore = 0;
      double exonScore = 0;

      Vector *newExons = Vector_new();

      // make an array containing cloned exons
      int k;
      for (k=0; k<Vector_getNumElement(model->features); k++) {
        SeqFeature *feature = Vector_getElementAt(model->features, k);

        if (feature->objectType != CLASS_DNADNAALIGNFEATURE) {
          Exon *origExon = (Exon *)feature;
          Exon *newExon = ExonUtils_cloneExon((Exon *)feature);
          // add in exon coverage scores from supporting features
          
          Vector *support = Exon_getAllSupportingFeatures(origExon);
          int m;
          for (m=0; m<Vector_getNumElement(support); m++) {
            BaseAlignFeature *daf = Vector_getElementAt(support, m);
            exonScore += BaseAlignFeature_getScore(daf);
          }
// NIY: Free supporting features vector??
          Exon_setStrand(newExon, strand);
          Vector_addElement(newExons, newExon);
        } else {
          Vector_addElement(newExons, feature);
        }
      }

      // trim the exons using the intron features to get the splice sites correct
// Note changed to only go up to num element - 1 because needs to be surrounded by exons and don't want to check in loop 
      for (k = 0; k<Vector_getNumElement(model->features)-1; k++) {
        SeqFeature *feature = Vector_getElementAt(model->features, k);
        if (feature->objectType == CLASS_DNADNAALIGNFEATURE) {
          DNAAlignFeature *intron = (DNAAlignFeature *)feature;
          if (DNAAlignFeature_getStrand(intron) != strand) {
            continue;
          }

          Exon *prevExon = Vector_getElementAt(newExons, k-1);
          Exon *nextExon = Vector_getElementAt(newExons, k+1);

          //next unless $new_exons[$k-1] && $new_exons[$k+1];
          // Do some sanity checking
          if (prevExon == NULL || nextExon == NULL ||
              prevExon->objectType != CLASS_EXON ||
              nextExon->objectType != CLASS_EXON) {
            fprintf(stderr, "Didn't get two exons surrounding intron\n");
            exit(EX_SOFTWARE);
          }

          Vector_addElement(introns, intron);

          // its an intron trim the exons accordingly
          Exon_setEnd(prevExon, DNAAlignFeature_getStart(intron));
          Exon_setStart(nextExon, DNAAlignFeature_getEnd(intron));
          if (Exon_getStart(prevExon) >= Exon_getEnd(prevExon)) { //$new_exons[$k-1]->start >=  $new_exons[$k-1]->end ) 
// Big jump - need to think about
//            next MODEL;
// Setting intron count to zero will cause the check after this loop to cause a 'continue' jump to the next model
            intronCount = 0;
            break;
          }

          intronCount++;
          intronScore += DNAAlignFeature_getScore(intron);
          if (DNAAlignFeature_getFlags(intron) & RSGINTRON_NONCANON) {
          //if (strstr(DNAAlignFeature_getHitSeqName(intron), "non canonical")) {
            nonConIntrons++;
          }

          // use the new intron feature code to store introns
          // provided we have not seen them before
          Intron *intronFeat;
          if (strand == 1) {
            intronFeat = Intron_new(prevExon, nextExon, DNAAlignFeature_getAnalysis(intron));
          } else {
            intronFeat = Intron_new(nextExon, prevExon, DNAAlignFeature_getAnalysis(intron));
          }
          IntronSupportingEvidence *ise = IntronSupportingEvidence_new();
          IntronSupportingEvidence_setValuesFromIntron(ise, intronFeat);
          IntronSupportingEvidence_setAnalysis(ise, DNAAlignFeature_getAnalysis(intron));
          IntronSupportingEvidence_setHitName(ise, DNAAlignFeature_getHitSeqName(intron));
          IntronSupportingEvidence_setScore(ise, DNAAlignFeature_getScore(intron));
          IntronSupportingEvidence_setScoreType(ise, "DEPTH");
          Intron_free(intronFeat);
    
          if (DNAAlignFeature_getFlags(intron) & RSGINTRON_NONCANON) {
          //if (strstr(DNAAlignFeature_getHitSeqName(intron), "non canonical")) {
            IntronSupportingEvidence_setIsSpliceCanonical(ise, 0);
          } else {
            IntronSupportingEvidence_setIsSpliceCanonical(ise, 1);
          }
          Vector_addElement(ises, ise);

          // Don't know why tested for ise defined push @ises, $ise if $ise;
        }
      }
      if (!intronCount) {
// Big jump - need to think about
        Vector_setFreeFunc(ises, IntronSupportingEvidence_freeImpl);
        Vector_free(ises);
        for (k=0;k<Vector_getNumElement(newExons); k++) {
          SeqFeature *sf = Vector_getElementAt(newExons,k);
          if (sf->objectType == CLASS_EXON) {
            Exon_freeImpl((Exon *)sf);
          }
        }
        Vector_free(newExons);
        Vector_free(introns);
        
        //next MODEL;
// NIY: Free stuff??
        //fprintf(stderr,"IntronCount continue\n");
        continue;
      }
      
      // trim padding from the start and end exons
      Exon *firstExon = Vector_getElementAt(newExons, 0);
      Exon_setStart(firstExon, Exon_getStart(firstExon) + 20);
      Exon *lastExon = Vector_getLastElement(newExons);
      Exon_setEnd(lastExon, Exon_getEnd(lastExon) - 20);

      // get rid of impossibly small exons
      // Note that here 'newExons' actually contains both 'introns' (DNAAlignFeatures) and Exons
      int tooShort = 0;
      for (k=0; k<Vector_getNumElement(newExons) && !tooShort; k++) {
        SeqFeature *sf = Vector_getElementAt(newExons, k);
        if (SeqFeature_getEnd(sf) - SeqFeature_getStart(sf) <= 0) {
// Big jump - need to think about
//          next MODEL;
          tooShort = 1;
        }
      }
      if (tooShort) {
        Vector_setFreeFunc(ises, IntronSupportingEvidence_freeImpl);
        Vector_free(ises);
        for (k=0;k<Vector_getNumElement(newExons); k++) {
          SeqFeature *sf = Vector_getElementAt(newExons,k);
          if (sf->objectType == CLASS_EXON) {
            Exon_freeImpl((Exon *)sf);
          }
        }
        Vector_free(newExons);
        Vector_free(introns);
// NIY: Free stuff??
        //fprintf(stderr,"TOOSHORT 1 continue\n");
        continue;
      }
      
      // trim away strings of Ns from the 1st and last exons
      // use same regex for 1st and last exon and reverse the
      // sequence accordingly depending on the strand
      Exon *newExon = Vector_getElementAt(newExons,0);
      char *exSeq = Exon_getSeqString(newExon);
      if (strand == -1) StrUtil_reverseString(exSeq, strlen(exSeq));

      int nN;
      for (nN=0; nN<Exon_getLength(newExon) && exSeq[nN] == 'N'; nN++) { }
          
      if (nN) {
        Exon_setStart(newExon, Exon_getStart(newExon) + nN);
      }

      newExon = Vector_getLastElement(newExons);
      exSeq = Exon_getSeqString(newExon);
      if (strand == 1) StrUtil_reverseString(exSeq, strlen(exSeq));

      for (nN=0; nN<Exon_getLength(newExon) && exSeq[nN] == 'N'; nN++) { }
          
      if (nN) {
        Exon_setEnd(newExon, Exon_getEnd(newExon) - nN);
      }
     
      // get rid of impossibly small exons again after the N trimming
      // Note that here 'newExons' actually contains both 'introns' (DNAAlignFeatures) and Exons
      for (k=0; k<Vector_getNumElement(newExons); k++) {
        SeqFeature *sf = Vector_getElementAt(newExons, k);
        if (SeqFeature_getEnd(sf) - SeqFeature_getStart(sf) <= 0) {
// Big jump - need to think about
//          next MODEL;
          tooShort = 1;
        }
      }
      if (tooShort) {
// NIY: Free stuff??
        //fprintf(stderr,"TOOSHORT 2 continue\n");
        continue;
      }

      // make it into a gene
      Vector *modifiedExons = Vector_new();
      for (k=0; k<Vector_getNumElement(newExons); k++) {
        SeqFeature *sf = Vector_getElementAt(newExons, k);
        if (sf->objectType != CLASS_DNADNAALIGNFEATURE) {
// Not sure if I really need to make these clones
          Vector_addElement(modifiedExons, ExonUtils_cloneExon((Exon *)sf));
          //Vector_addElement(modifiedExons, sf);
   
          Object_decRefCount(sf);
          Exon_free(sf);
        }
      }
      Vector_free(newExons);
 
      if ( strand == 1 ) {
        Vector_sort(modifiedExons, SeqFeature_startCompFunc);
      } else {
        Vector_sort(modifiedExons, SeqFeature_reverseStartCompFunc);
      }

      // make it into a gene
//      my $t =  new Bio::EnsEMBL::Transcript(-EXONS => \@modified_exons);
      Transcript *t = Transcript_new();
      for (k=0; k<Vector_getNumElement(modifiedExons); k++) {
        Exon *exon = Vector_getElementAt(modifiedExons, k);

        Transcript_addExon(t, exon, 0);// k+1);
      }
      Vector_free(modifiedExons);
 
      //fprintf(stderr,"Adding ISEs to transcript\n");
      for (k=0; k<Vector_getNumElement(ises); k++) {
        IntronSupportingEvidence *ise = Vector_getElementAt(ises, k);
        Transcript_addIntronSupportingEvidence(t, ise);
      }
// NIY: Free ises vector??
      Vector_free(ises);

//# check for dna
//#      my $check = $t->seq->seq ;
//#      my $Ns =  $check =~  s/N//g;
//#      if( length($t->seq->seq) == $Ns ){
//#        $self->throw("There does not appear to be ay DNA in the database, transcript seq is all N's\n");
//#      }
      // add a translation 
      // NOTE: Extra cloning my $tran = compute_translation(clone_Transcript($t));
      Transcript *tran = TranslationUtils_computeTranslation(t);

      // stop spam coming from the Exon module
      Transcript_setDbID(tran, 0);
      // keep track of the scores for this transcript
      Transcript_setAnalysis(tran, RefineSolexaGenes_getAnalysis(rsg));
      Transcript_setVersion(tran, 1);

      // favor longer cds by adding doubling the score for coding exons
      // only use exons that are completely coding otherwise you also 
      // end up adding in score which is really UTR for long terminal exons
      // that have a bit of coding in them
      double codingBonus = 0;
      int codingExons = 0;
      if ( Transcript_getTranslation(tran) != NULL ) {
        for (k=0; k<Transcript_getExonCount(tran); k++) {
          Exon *ce = Transcript_getExonAt(tran, k);
          if (! (Exon_getPhase(ce) == -1 || Exon_getEndPhase(ce) == -1)) {
            Vector *support = Exon_getAllSupportingFeatures(ce);
            if (Vector_getNumElement(support) > 0) {
              DNAAlignFeature *daf = Vector_getElementAt(support, 0);
  
              codingBonus += DNAAlignFeature_getScore(daf);
            } else {
              if (verbosity > 0) fprintf(stderr, "ERROR ERROR ERROR: Have an exon with no support: Exon %p\n", ce);
            }
            codingExons++;
// NIY: Free support?
          }
        }
      }
      
      // remove any supporting features from the exons
      for (k=0; k<Transcript_getExonCount(tran); k++) {
        Exon *e = Transcript_getExonAt(tran, k);

        Exon_flushSupportingFeatures(e);
      }

// There's a bunch of hackery here adding random stuff to the transcript

      // print "Coding Bonus of $coding_bonus from $coding_exons completely coding exons \n";
      Transcript_setScore(tran,(((int)( intronScore + exonScore )) / 10.0 ) + codingBonus);

      TranscriptExtraData *ted = TranscriptExtraData_new();

      // store the introns along with the transcript so we can use them later for UTR trimming 
      // Moved this line down a bit to be near other hackery
      ted->introns = introns; // Note introns are infact DNAAlignFeatures here !!!!

      // print "Final score = $intron_score + int( $exon_score / 100 ) + $coding_bonus = " . $tran->{'_score'} ;
// Doesn't seem to be used      $tran->{'_depth'} =  ( $intron_score + $exon_score );
      ted->depth = intronScore + exonScore;
// Only used in this method
      ted->nNCIntrons =  nonConIntrons;
// Used in another method
      ted->exonUse = StrUtil_copyString(&ted->exonUse, exonUse, 0);
// Only used in this method
      ted->intronCount = intronCount;

      Transcript_setExtraData(tran, ted);
      Vector_addElement(trans, tran);

      // we want X number of models
      if (RefineSolexaGenes_getBestScoreType(rsg) && RefineSolexaGenes_getMaxNum(rsg)) {
        if (Vector_getNumElement(trans) > RefineSolexaGenes_getMaxNum(rsg)) {
//NIY: Possibly some tidying needed here!
          break;
    //      last MODEL if scalar(@trans)  >= ( $self->MAX_NUM +1 )  ;
        }
      }
    Translation *tmp_tln = Transcript_getTranslation(tran);
    Vector *exons = Transcript_getAllExons(tran);
    Exon *startExon = Translation_getStartExon(tmp_tln);
    Exon *endExon = Translation_getEndExon(tmp_tln);
    int m;
    int startRank;
    int endRank;
    for (m = 0; m < Vector_getNumElement(exons); m++) {
      if (startExon == Vector_getElementAt(exons, m)) {
        startRank = m;
      }
      if (endExon == Vector_getElementAt(exons, m)) {
        endRank = m;
        break;
      }
    }
    fprintf(stderr, "TIBO S %d %ld %ld E %d %ld %ld\n", startRank, Exon_getStart(startExon), Exon_getStart(startExon), endRank, Exon_getStart(endExon), Exon_getStart(endExon));
    }

    Vector_free(modelsByScore);

    // re-sort the transcripts to take account of the revised scores
    Vector_sort(trans, SeqFeature_reverseScoreCompFunc);
    //@trans = sort { $b->{'_score'} <=> $a->{'_score'} } @trans;

    if (Vector_getNumElement(trans) && !cluster->finalModels) cluster->finalModels = Vector_new();
// Doesn't seem to be used (only set)    Transcript *best;
    //fprintf(stderr,"Have %d transcripts\n", Vector_getNumElement(trans));
    for (j=0; j<Vector_getNumElement(trans); j++) {
      Transcript *tran = Vector_getElementAt(trans, j);

      TranscriptExtraData *ted = Transcript_getExtraData(tran);

      Gene *newGene = TranscriptUtils_convertToGene(tran, Gene_getAnalysis(gene), NULL);

      version++;
      Gene_setBiotype(newGene, RefineSolexaGenes_getOtherIsoformsType(rsg));

      if (version == 1) {
        Gene_setBiotype(newGene, RefineSolexaGenes_getBestScoreType(rsg));
        //fprintf(stderr,"Best model set to gene at %ld %ld %d\n",Gene_getStart(newGene), Gene_getEnd(newGene), Gene_getStrand(newGene)); 
//        best = tran;
      }


      char stableId[2048];
      sprintf(stableId, "%s-v%d.%d-%d-%d-%d-NC-%d-%d",
              Gene_getStableId(gene),
              clusterCount,
              version,
              (int)(Transcript_getScore(tran)),
              ted->depth,
              ted->intronCount,
              ted->nNCIntrons, 
              Transcript_getStrand(tran));
      Gene_setStableId(newGene, stableId);
      Vector_addElement(cluster->finalModels, newGene);
    }
    Vector_free(trans);
  }

  //fprintf(logfp, "Done gene cycle\n");
  return modelClusters;
}


// Hacked version of method from TranscriptUtils
Gene *TranscriptUtils_convertToGene(Transcript *t, Analysis *analysis, char *biotype) {
  Gene *g;

// Ignore for now
//  fullyLoadTranscript(t);

  if (analysis == NULL)  analysis = Transcript_getAnalysis(t);
  if (biotype != NULL) Transcript_setBiotype(t, biotype);
 
  g = Gene_new();
  if (Transcript_getBiotype(t)) Gene_setBiotype   (g, Transcript_getBiotype(t));
  Gene_addTranscript(g, t);
  Gene_setAnalysis  (g, analysis);
// Odd!!!
  Gene_setDbID      (g, Transcript_getDbID(t));

  Gene_setStrand(g, Transcript_getStrand(t));
  Gene_setStart(g, Transcript_getStart(t));
  Gene_setEnd(g, Transcript_getEnd(t));

/*
  fprintf(stderr,"converted transcript with coords %ld-%ld:%d to gene with coords %ld-%ld:%d\n",
          Transcript_getStart(t), Transcript_getEnd(t), Transcript_getStrand(t),
          Gene_getStart(g), Gene_getEnd(g), Gene_getStrand(g));
*/
  
/* Wierd checks
    throw("there are no transcripts for this gene\n") if scalar(@{$g->get_all_Transcripts}) == 0 ;
    for my $tr ( @{$g->get_all_Transcripts} ) {
      throw("there are no exons  for this transcript \n") if scalar(@{$tr->get_all_Exons}) == 0 ;
    }
    throw("Problems with ".Gene_info($gene)." undef coords") if(!$gene->start || !$gene->end);
*/
  return g;
}


//We return 0 if the gene is OK, 1 if the gene failed a test
int RefineSolexaGenes_pruneUTR(RefineSolexaGenes *rsg, Gene *gene) {
  int verbosity = RefineSolexaGenes_getVerbosity(rsg);
  //fprintf(stderr, "pruneUTR called\n");
  if (!RefineSolexaGenes_trimUTR(rsg)) {
    return 0;
  } 

  Transcript *transcript = Gene_getTranscriptAt(gene, 0);
  if (!Transcript_getTranslation(transcript)) {
    if (verbosity > 0) fprintf(stderr,"No translation so returning gene unchanged\n");
    return 0;
  }
 
  // fetch introns 
  TranscriptExtraData *ted = Transcript_getExtraData(transcript);

  Vector *introns = ted->introns;

  // otherwise trim the UTR according to the values set out in the config
  StringHash *intronHash = StringHash_new(STRINGHASH_SMALL);

  int i;
  char key[1024];
  for (i=0; i<Vector_getNumElement(introns); i++) {
    SeqFeature *sf = Vector_getElementAt(introns, i);
    
    sprintf(key,"%ld:%ld:%d", SeqFeature_getStart(sf), SeqFeature_getEnd(sf), SeqFeature_getStrand(sf));
    if (!StringHash_contains(intronHash, key)) {
      StringHash_add(intronHash, key, sf);
    } else {
      if (verbosity > 0) fprintf(stderr,"Error: got an intron (Dna) feature twice with key %s\n", key);
    }
  }
  
  Vector *newFiveP = Vector_new();
  Vector *newThreeP = Vector_new();
  Vector *features = Vector_new();

  Vector *exons = Vector_copy(Transcript_getAllExons(transcript));
  Vector_sort(exons, SeqFeature_startCompFunc);
  
  // put everything into the features array
  Vector_append(features, exons);

  for (i=0; i<Vector_getNumElement(exons)-1; i++) {
    Exon *exon = Vector_getElementAt(exons, i); 
    Exon *exonP1 = Vector_getElementAt(exons, i+1); 
    sprintf(key, "%ld:%ld:%d", Exon_getEnd(exon), Exon_getStart(exonP1), Exon_getStrand(exon));

    if (StringHash_contains(intronHash, key)) {
      SeqFeature *intron = StringHash_getValue(intronHash, key);
      Vector_addElement(features, intron);
    } else {
      if (verbosity > 0) fprintf(stderr,"Error: couldn't find an intron (Dna) feature with key %s\n", key);
    }
  }
  StringHash_free(intronHash, NULL);

  Vector_sort(features, SeqFeature_startCompFunc);
  // so now we should have an array of alternating introns and exons
  if (verbosity > 0) fprintf(logfp, "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\nTrimming UTR\nTranscript %s %ld %ld %d %d\n\n",
         Transcript_getSeqRegionName(transcript), Transcript_getStart(transcript), Transcript_getEnd(transcript), Transcript_getStrand(transcript),
         Transcript_getExonCount(transcript));

  if (Vector_getNumElement(features) != (Vector_getNumElement(exons)*2)-1) {
    fprintf(stderr,"Something is wrong we have wrong number of introns compared to exons"
                   " - have %d introns, %d exons and %d features in combined vector\n", 
            Vector_getNumElement(introns), Vector_getNumElement(exons), Vector_getNumElement(features));
    exit(EX_SOFTWARE);
  }
  double averageIntron = 0;
  int intronCount = 0;
  Translation *translation = Transcript_getTranslation(transcript);
  // leave single exon genes alone for now
  if (Vector_getNumElement(features) == 1 || Translation_getStartExon(translation) == Translation_getEndExon(translation)) { // was: scalar(@{$transcript->get_all_translateable_Exons}) == 1 ) 
    // lets strip the UTR
// NIY: Probably a leak here
    //fprintf(stderr, "Single exon transcript\n");
    Vector *translateableExons = Transcript_getAllTranslateableExons(transcript);

    Vector *newExons = Vector_new();
    Vector_addElement(newExons, ExonUtils_cloneExon(Vector_getElementAt(translateableExons, 0)));

    Transcript_freeAdjustedTranslateableExons(transcript, translateableExons);
    Vector_free(translateableExons);

    Transcript *trimmedTran = RefineSolexaGenes_modifyTranscript(rsg, transcript, newExons);

    Vector_free(newExons);

    // The naughty bit!
// NIY: Definitely a leak here
    Vector_free(gene->transcripts);

    gene->transcripts = Vector_new();
    Gene_addTranscript(gene, trimmedTran);

    Vector_free(features);
    Vector_free(exons);
    return 0;
  }

  // first calculate the average
  for (i=0; i<Vector_getNumElement(features); i++) {
    SeqFeature *f = Vector_getElementAt(features, i);
    if (f->objectType == CLASS_DNADNAALIGNFEATURE) {
      averageIntron += SeqFeature_getScore(f);
      intronCount++;
    }
  }
  averageIntron /= intronCount;

  for (i=0; i<Vector_getNumElement(features); i++) {
    SeqFeature *f = Vector_getElementAt(features, i);
    if (f->objectType == CLASS_DNADNAALIGNFEATURE) {
      if (averageIntron && (SeqFeature_getScore(f) / averageIntron) * 100.0 <= RefineSolexaGenes_getRejectIntronCutoff(rsg)) {
        if (verbosity > 0) fprintf(logfp, " Potentially bad ");
      }
    }
  }

  if (intronCount != Vector_getNumElement(exons)-1) {
    fprintf(stderr, "Error: Something is wrong we are missing introns %d exons  and %d introns\n", Vector_getNumElement(exons), intronCount);
    exit(EX_SOFTWARE);
  }
  //fprintf(logfp, "Average intron depth = %lf\n", averageIntron);

  Vector_free(exons);
  
  Vector *fiveP  = Vector_new();
  Vector *threeP = Vector_new();
  int coding = 0;

  // need to account for strand
  if (Transcript_getStrand(transcript) == -1) {
    Vector_sort(features, SeqFeature_reverseStartCompFunc);
  }

  for (i=0; i<=Vector_getNumElement(features); i+=2) {
    Exon *e = Vector_getElementAt(features, i);
    if (e->objectType != CLASS_EXON) {
      fprintf(stderr, "Got a DNA align feature where I should have an exon\n");
      exit(EX_SOFTWARE);
    }

// Perl did    if ( $e->coding_region_start($transcript) ) 
    if (Translation_getStartExon(translation) == e) {
      // first coding exon
      int j;
      for (j=0 ; j<=i; j++) {
        Vector_addElement(fiveP, Vector_getElementAt(features, j));
      }
      break;
    }
  }

  for (i=Vector_getNumElement(features)-1; i>=0; i-=2) {
    Exon *e = Vector_getElementAt(features, i);
    if (e->objectType != CLASS_EXON) {
      fprintf(stderr, "Got a DNA align feature where I should have an exon\n");
      exit(EX_SOFTWARE);
    }
// Perl was   if ( $e->coding_region_end($transcript) ) 
    if (Translation_getEndExon(translation) == e) {
      // first (last) coding exon
      int j;
      for (j=i ; j<Vector_getNumElement(features); j++) {
        Vector_addElement(threeP, Vector_getElementAt(features, j));
      }
      break;
    }
  }
  Vector_free(features);
  
  // want to start at last coding exon and work outwards so....
  Vector_reverse(fiveP);

  // now we should be good
  if (verbosity > 0) fprintf(logfp, "FIVE P\n");

// NIY Not sure if need to make a copy here (perl did)
  Vector *newExons = Transcript_getAllTranslateableExons(transcript);
  //fprintf(stderr, "num element in newExons before removing start cds and end cds exons = %d\n", Vector_getNumElement(newExons));

  Exon *fivePCDS  = Vector_removeElementAt(newExons, 0);
  Exon *threePCDS = Vector_removeElementAt(newExons, Vector_getNumElement(newExons)-1);

  //fprintf(stderr, "num element in newExons after removing start cds and end cds exons = %d\n", Vector_getNumElement(newExons));
//Perl was:
//  my $fivep_cds = shift(@new_exons);
//  my $threep_cds = pop(@new_exons);

  long fivePLen;
  long threePLen;
  int fivePc  = 0;
  int threePc = 0;
  long nmd;
  
  // FIVE PRIME RULES
  
// FIVEP: 
  for (i = 0; i < Vector_getNumElement(fiveP); i++) {
    SeqFeature *f = Vector_getElementAt(fiveP, i);

    if (verbosity > 0) fprintf(logfp, "%ld\t%ld\t%s",SeqFeature_getStart(f), SeqFeature_getEnd(f), Class_findByType(f->objectType)->name);

    if (i == 0) {
      if (f->objectType != CLASS_EXON) {
        fprintf(stderr, "First feature is not an exon\n");
        exit(EX_SOFTWARE);
      }

      // UTR starts in this exon - how long is it?
      long cdsStart = Exon_getCodingRegionStart((Exon *)f, transcript);
      if (Transcript_getStrand(transcript) == -1) {
        cdsStart = Exon_getCodingRegionEnd((Exon *)f, transcript);
      }

      if (cdsStart == POS_UNDEF) {
        fprintf(stderr, "First coding exon has no CDS\n");
        exit(EX_SOFTWARE);
      }

      if (verbosity > 0) fprintf(logfp, "\tCDS START %ld\t", cdsStart);

      if (Transcript_getStrand(transcript) == 1) {
        fivePLen = cdsStart - SeqFeature_getStart(f) + 1;
      } else {
        fivePLen = SeqFeature_getEnd(f) - cdsStart + 1;
      }

      // is the coding exon too long
      if (fivePLen > RefineSolexaGenes_getMax5PrimeLength(rsg)) {
        // replace it with the cds
// Not sure there's anything in newFiveP at this point
        Vector_removeAll(newFiveP);

        Vector_addElement(newFiveP, fivePCDS);
        if (verbosity > 0) fprintf(logfp, " 5p too long %ld\n", fivePLen);
        //last FIVEP;
        break;
      }
      Vector_addElement(newFiveP, f);
      fivePc++;

    } else {

      if (f->objectType == CLASS_EXON) {
        fivePc++;
        fivePLen += SeqFeature_getLength(f);

        // does it make the UTR too long?
        if (fivePLen > RefineSolexaGenes_getMax5PrimeLength(rsg)) {
          // dont add it
          if (verbosity > 0) fprintf(logfp, " 5p too long %ld\n", fivePLen);
          //last FIVEP;
          break;
        }
        // is it too many exons?
        if (fivePc > RefineSolexaGenes_getMax5PrimeExons(rsg)) {
          // dont add it
          if (verbosity > 0) fprintf(logfp, " too many 5p  %d cut them all as we are not sure\n", fivePc);

          Vector_removeAll(newFiveP);

          Vector_addElement(newFiveP, fivePCDS);
          //last FIVEP;
          break;
        }
        Vector_addElement(newFiveP, f);
      }
    }
    // Does the intron score well enough to include the exon
    // apply rules and add successful exons into the mix
    if (f->objectType == CLASS_DNADNAALIGNFEATURE) {
      if (verbosity > 0) fprintf(logfp, "\t%lf",SeqFeature_getScore(f));

      if (averageIntron && (SeqFeature_getScore(f) / averageIntron) * 100.0 <= RefineSolexaGenes_getRejectIntronCutoff(rsg)) {
        if (verbosity > 0) fprintf(logfp, " Rejecting on score cutoff %f vs %lf\n", SeqFeature_getScore(f), averageIntron);
        // dont add any more 
        //last FIVEP;
        break;
      }
    }
    if (verbosity > 0) fprintf(logfp, "\n");
  }
  Vector_free(fiveP);
  
  if (verbosity > 0) fprintf(logfp, "\n");
  // three P
  if (verbosity > 0) fprintf(logfp, "THREE P\n");

// THREEP:   
  long new_exonEnd;
  long end;
  long start;
  for (i = 0; i < Vector_getNumElement(threeP); i++) {
    SeqFeature *f = Vector_getElementAt(threeP, i);

    if (verbosity > 0) fprintf(logfp, "%ld\t%ld\t%s",SeqFeature_getStart(f), SeqFeature_getEnd(f), Class_findByType(f->objectType)->name);

    if (i == 0) {
      if (f->objectType != CLASS_EXON) {
        fprintf(stderr, "First feature is not an exon\n");
        exit(EX_SOFTWARE);
      }

      // UTR starts in this exon - how long is it?
      long cdsEnd = Exon_getCodingRegionEnd((Exon *)f, transcript);
      if (Transcript_getStrand(transcript) == -1) {
        cdsEnd = Exon_getCodingRegionStart((Exon *)f, transcript);
      }

      if (cdsEnd == POS_UNDEF) {
        fprintf(stderr, "Last coding exon has no CDS\n");
        exit(EX_SOFTWARE);
      }

      if (verbosity > 0) fprintf(logfp, "\tCDS END %ld\t", cdsEnd);
      if (SeqFeature_getStrand(f) == 1) {
        start = cdsEnd;
        end = SeqFeature_getEnd(f);
      }
      else {
        start = SeqFeature_getStart(f);
        end = cdsEnd;
      }
      new_exonEnd = RefineSolexaGenes_has_polyA_signal(rsg, (Exon *)f, start, end);

      if (Transcript_getStrand(transcript) == -1) {
        threePLen = cdsEnd - SeqFeature_getStart(f) + 1;
      } else {
        threePLen = SeqFeature_getEnd(f) - cdsEnd + 1;
      }

      // is the coding exon too long
      if (threePLen > RefineSolexaGenes_getMax3PrimeLength(rsg)) {
          if (new_exonEnd) {
            if (Transcript_getStrand(transcript) == -1) {
              threePLen = cdsEnd - new_exonEnd + 1;
            } else {
              threePLen = new_exonEnd - cdsEnd + 1;
            }
            if (threePLen < RefineSolexaGenes_getMax3PrimeLength(rsg)) {
              if (SeqFeature_getStrand(f) == 1) {
                SeqFeature_setEnd(f, new_exonEnd);
              }
              else {
                SeqFeature_setStart(f, new_exonEnd);
              }
              fprintf(logfp, "Cutting current at %ld to avoid - rejected 3' length\n", new_exonEnd);
              Vector_addElement(newThreeP, f);
              break;
            }
          }
        // replace it with the cds
// Not sure there's anything in newThreeP at this point
        Vector_removeAll(newThreeP);

        Vector_addElement(newThreeP, threePCDS);
        if (verbosity > 0) fprintf(logfp, " 3p too long %ld\n", threePLen);
        //last THREEP;
        break;
      }
      Vector_addElement(newThreeP, f);
      nmd = threePLen ;
      threePc++;

    } else {

      if (f->objectType == CLASS_EXON) {
        // does it break the NMD rule?
        if (nmd > 55) {
          if (verbosity > 0) fprintf(logfp, " splice is after %ld bp from stop codon - rejected on NMD rule of maximum 55 bp\n", nmd);
          if (new_exonEnd) {
            fprintf(logfp, "Cutting previous at %ld to avoid - rejected on NMD rule of maximum 55 bp %ld\n", new_exonEnd, nmd);
            SeqFeature *feat = Vector_getElementAt(newThreeP, 0);
            if (SeqFeature_getStrand(feat) == 1) {
              SeqFeature_setEnd(feat, new_exonEnd);
            }
            else {
              SeqFeature_setStart(feat, new_exonEnd);
            }
          }
          else {
            Vector_removeAll(newThreeP);

            Vector_addElement(newThreeP, threePCDS);
          }
          //last THREEP;
          break;
        }
        if (SeqFeature_getStrand(f) == 1) {
          start = SeqFeature_getStart(f);
          end = SeqFeature_getEnd(f);
        }
        else {
          start = SeqFeature_getStart(f);
          end = SeqFeature_getEnd(f);
        }
        // is it too many exons?
        if (threePc > RefineSolexaGenes_getMax3PrimeExons(rsg)) {
          // dont add it
          if (verbosity > 0) fprintf(logfp, " too many 3p  %d cut them all as we are not sure\n", threePc);
          if (new_exonEnd) {
            SeqFeature *feat = Vector_getElementAt(newThreeP, 0);
            if (SeqFeature_getStart(f) < new_exonEnd) {
              if (SeqFeature_getStrand(feat) == 1) {
                SeqFeature_setEnd(feat, new_exonEnd);
              }
              else {
                SeqFeature_setStart(feat, new_exonEnd);
              }
              fprintf(logfp, "Cutting current at %ld to avoid - rejected 3' number of exon\n", new_exonEnd);
              break;
            }
          }

          Vector_removeAll(newThreeP);

          Vector_addElement(newThreeP, threePCDS);
          //last THREEP;
          break;
        }
        new_exonEnd = RefineSolexaGenes_has_polyA_signal(rsg, (Exon *)f, start, end);
        threePc++;
        threePLen += SeqFeature_getLength(f);

        // does it make the UTR too long?
        if (threePLen > RefineSolexaGenes_getMax3PrimeLength(rsg)) {
          // dont add it
          if (verbosity > 0) fprintf(logfp, " 3p too long %ld\n", threePLen);
          if (new_exonEnd) {
            threePLen -= SeqFeature_getLength(f);
            if (Transcript_getStrand(transcript) == -1) {
              threePLen += SeqFeature_getEnd(f) - new_exonEnd + 1;
            } else {
              threePLen = new_exonEnd - SeqFeature_getStart(f) + 1;
            }
            if (threePLen < RefineSolexaGenes_getMax3PrimeLength(rsg)) {
              fprintf(logfp, "Cutting current at %ld to avoid - rejected 3' length\n", new_exonEnd);
              if (SeqFeature_getStrand(f) == 1) {
                SeqFeature_setEnd(f, new_exonEnd);
              }
              else {
                SeqFeature_setStart(f, new_exonEnd);
              }
            }
            Vector_addElement(newThreeP, f);
            break;
          }
          //last THREEP;
          break;
        }
        Vector_addElement(newThreeP, f);
      }
    }

    // Does the intron score well enough to include the exon
    // apply rules and add successful exons into the mix
    if (f->objectType == CLASS_DNADNAALIGNFEATURE) {
      if (verbosity > 0) fprintf(logfp, "\t%lf",SeqFeature_getScore(f));

      if (averageIntron && (SeqFeature_getScore(f) / averageIntron) * 100.0 <= RefineSolexaGenes_getRejectIntronCutoff(rsg)) {
        if (verbosity > 0) fprintf(logfp, " Rejecting on score cutoff %f vs %lf\n", SeqFeature_getScore(f), averageIntron);
        // dont add any more
        //last THREEP;
        break;
      }
    }
    if (verbosity > 0) fprintf(logfp, "\n");
  }
  Vector_free(threeP);
  
  Vector_append(newExons, newFiveP);
  Vector_append(newExons, newThreeP);

  Vector_free(newFiveP);
  Vector_free(newThreeP);

  if (verbosity > 0) fprintf(logfp, " New transript has %d exons\n", Vector_getNumElement(newExons));

  Vector *clones = Vector_new();

  for (i=0; i<Vector_getNumElement(newExons); i++) {
    Exon *e = Vector_getElementAt(newExons, i);
    if (e->objectType != CLASS_EXON) {
      fprintf(stderr, "Feature is not an exon (start %ld end %ld %s (%p)\n",
              SeqFeature_getStart(e), SeqFeature_getEnd(e), Class_findByType(e->objectType)->name, e);
      exit(EX_SOFTWARE);
    }
    Vector_addElement(clones, ExonUtils_cloneExon(e));
  }

  if (Transcript_getStrand(transcript) == 1) {
    Vector_sort(clones, SeqFeature_startCompFunc);
  } else {
    Vector_sort(clones, SeqFeature_reverseStartCompFunc);
  }

  Transcript *trimmedTran = RefineSolexaGenes_modifyTranscript(rsg, transcript, clones);

  Vector_free(clones);

  // The naughty bit!
// NIY: Definitely a leak here
  Vector_free(gene->transcripts);

  gene->transcripts = Vector_new();
  Gene_addTranscript(gene, trimmedTran);
  // filter single exon genes that may have been made through UTR trimming
//  if (Gene_getExonCount(gene) == 1) {
//    if (RefineSolexaGenes_getSingleExonModelType(rsg)) {
//      Exon *exon = Transcript_getExonAt(Gene_getTranscriptAt(gene, 0), 0);
//      if (Exon_getLength(exon) < RefineSolexaGenes_getMinSingleExonLength(rsg)) {
//        if (verbosity > 0) fprintf(stderr, "Skipping single exon gene after UTR pruning - length less than MinSingleExonLength\n");
//        //We don't want to store the gene so we return 1
//        return 1;
//      }
//      Gene_setBiotype(gene, RefineSolexaGenes_getSingleExonModelType(rsg));
//    } else {
//      if (verbosity > 0) fprintf(stderr, "Skipping single exon gene after UTR pruning - no SingleExonModelType defined\n");
//      //We don't want to store the gene so we return 1
//      return 1;
//    }
//  }

// NIY: Free up everything!

  return 0;
}


// NOTE: Exons MUST BE sorted in correct order for adding to transcript
Transcript *RefineSolexaGenes_modifyTranscript(RefineSolexaGenes *rsg, Transcript *tran, Vector *exons) {
  long cdsStart;
  long cdsEnd;
  int verbosity = RefineSolexaGenes_getVerbosity(rsg);

  if (Transcript_getStrand(tran) == 1) {
    cdsStart = Transcript_getCodingRegionStart(tran);
    cdsEnd   = Transcript_getCodingRegionEnd(tran);
  } else {
    cdsEnd   = Transcript_getCodingRegionStart(tran);
    cdsStart = Transcript_getCodingRegionEnd(tran);
  }
  
  if (verbosity > 0) fprintf(logfp, "CDS START END %ld  %ld \n", cdsStart,  cdsEnd);
  Translation *tln = Transcript_getTranslation(tran); 
  if (verbosity > 0) fprintf(logfp, "PHASE %d %d\n", Translation_getStart(tln), Translation_getEnd(tln));

// t created here and again below - not sure if need both
  Transcript *t = Transcript_new();
  int i;
  for (i=0; i<Vector_getNumElement(exons); i++) {
    Exon *exon = Vector_getElementAt(exons, i);
    Transcript_addExon(t, exon, 0);//i+1);
  }

  Exon *se = NULL;
  Exon *ee = NULL;
  for (i=0; i<Transcript_getExonCount(t); i++) {
    Exon *e = Transcript_getExonAt(t, i);
    //fprintf(stderr,"looking for ee and se have exon %ld %ld cds start %ld cds end %ld\n", 
    //        Exon_getStart(e), Exon_getEnd(e), cdsStart, cdsEnd);
    if (Exon_getStart(e) <= cdsEnd && Exon_getEnd(e) >= cdsEnd) {
     // fprintf(stderr,"found ee\n");
      ee = e;
    }
    if (Exon_getStart(e) <= cdsStart && Exon_getEnd(e) >= cdsStart) {
      //fprintf(stderr,"found se\n");
      se = e;
    }
  }

  if (ee == NULL || se == NULL) {
    fprintf(stderr,"NULL for ee or se\n");
    exit(EX_SOFTWARE);
  }
  long ts;
  long te;
  if ( Transcript_getStrand(tran) == -1 ) {
    //fprintf(stderr,"reverse strand in modifyTranscript\n");
    ts =  Exon_getEnd(se) - cdsStart + 1;
    te =  Exon_getEnd(ee) - cdsEnd   + 1;
  } else {
    //fprintf(stderr,"forward strand in modifyTranscript\n");
    ts =   cdsStart - Exon_getStart(se) + 1;
    te =   cdsEnd   - Exon_getStart(ee) + 1;
  }

// Why created again?????
/*
  Transcript *t = Transcript_new();
  Transcript_setExons(exons);
*/

  // transfer the intron supporting evidence
  // except for where we have trimmed the intron
  Vector *intronSupport = Transcript_getAllIntronSupportingEvidence(tran);
  if (intronSupport!=NULL) {
    for (i=0; i<Vector_getNumElement(intronSupport); i++) {
      IntronSupportingEvidence *ise = Vector_getElementAt(intronSupport, i);
      if (IntronSupportingEvidence_getSeqRegionStart(ise) > Transcript_getStart(t) &&  
          IntronSupportingEvidence_getSeqRegionEnd(ise) < Transcript_getEnd(t)) {
        Transcript_addIntronSupportingEvidence(t, ise);
        //Object_decRefCount(ise);
      }
    }
  }
  //#  my $start_phase = $se->phase;
  Translation *translation = Translation_new();
  Translation_setStartExon(translation, se); 
  Translation_setEndExon(translation, ee); 
  Translation_setStart(translation, ts); 
  Translation_setEnd(translation, te); 

  if (verbosity > 0) fprintf(logfp,"S-E %ld %ld \n",ts, te); //#START PHASE $start_phase\n";
  if (verbosity > 0) fprintf(logfp,"GS %ld %ld\n", Translation_getGenomicStart(translation), Translation_getGenomicEnd(translation));
  Transcript_setTranslation(t, translation);
  //# calculate_exon_phases($t,$start_phase);

  char *tranTranslationSeq = Transcript_translate(tran);
  char *tTranslationSeq = Transcript_translate(t);
  if (strcmp(tranTranslationSeq, tTranslationSeq)) {
    fprintf(stderr, "Translations do not match: Before %s\nAfter  %s\n", Transcript_translate(tran), Transcript_translate(t));
    exit(EX_SOFTWARE);
  }
// NIY: free translation strings??
  free(tranTranslationSeq);
  free(tTranslationSeq);
  if (tran->extraData) TranscriptExtraData_free(tran->extraData);
  Transcript_free(tran);
  return t;
}

bam_hdr_t *RefineSolexaGenes_getAlignmentBamHeader(RefineSolexaGenes *rsg) {
  if (rsg->alignment_header == NULL) {
    char *bamFile = RefineSolexaGenes_getAlignmentBamFile(rsg);
    rsg->alignment_sam = hts_open(bamFile, "rb");
    if (rsg->alignment_sam == NULL) {
      fprintf(stderr, "Bam file %s not found\n", bamFile);
      exit(EX_NOINPUT);
    }
    fprintf(stderr,"Opened bam file %s\n", bamFile);

    rsg->alignment_idx = sam_index_load(rsg->alignment_sam, bamFile); // load BAM index
    if (rsg->alignment_idx == 0) {
      fprintf(stderr, "BAM index file is not available.\n");
      exit(EX_NOINPUT);
    }
    fprintf(stderr,"Opened bam index for %s\n", bamFile);

    rsg->alignment_header = bam_hdr_init();
    rsg->alignment_header = bam_hdr_read(rsg->alignment_sam->fp.bgzf);
  }
  return rsg->alignment_header;
}

hts_idx_t *RefineSolexaGenes_getAlignmentBamIndex(RefineSolexaGenes *rsg) {
  return rsg->alignment_idx;
}

void RefineSolexaGenes_cleanPolyAcheck(RefineSolexaGenes *rsg) {
  if (rsg->alignment_header != NULL) {
    hts_idx_destroy(rsg->alignment_idx);
    bam_hdr_destroy(rsg->alignment_header);
    hts_close(rsg->alignment_sam);
  }
}

long RefineSolexaGenes_has_polyA_signal(RefineSolexaGenes *rsg, Exon *exon, long start, long end) {
    fprintf(stderr, "Working on %ld %ld\n", start, end);
  Slice *slice = Exon_getSlice(exon);
  CachingSequenceAdaptor *cachingSeqAdaptor = DBAdaptor_getCachingSequenceAdaptor(Slice_getAdaptor(slice)->dba);
  int strand = Exon_getStrand(exon);
  char *seq = CachingSequenceAdaptor_fetchBySliceStartEndStrand(cachingSeqAdaptor, 
                                                                          slice,
                                                                          start,
                                                                          end,
                                                                          strand);
  char region[2048];
  char region_name[1024];
  int ref;
  int begRange;
  int endRange;
  long cdsend = 0;
  long pos = 0;
  long lastpos = 0;
  char *res = strstr(seq, "AATAAA");
  if (RefineSolexaGenes_getUcscNaming(rsg) == 0) {
    sprintf(region_name,"%s", Slice_getSeqRegionName(slice));
  } else {
    sprintf(region_name,"chr%s", Slice_getSeqRegionName(slice));
  }
  bam_hdr_t *header = RefineSolexaGenes_getAlignmentBamHeader(rsg);
  ref = bam_name2id(header, region_name);
  if (ref < 0) {
    fprintf(stderr, "Invalid region %s\n", region_name);
    exit(EX_DATAERR);
  }
  bam1_t *read = bam_init1();
  while (res != NULL) {
    pos = (long)(res-seq)+start;
    fprintf( stderr, "Found AATAAA at position %ld\n%s\n", pos, res);
    sprintf(region,"%s:%ld-%ld", region_name, pos,
                               pos+40);

    if (hts_parse_reg(region, &begRange, &endRange) == NULL) {
      fprintf(stderr, "Could not parse %s\n", region);
      exit(EX_DATAERR);
    }
    fprintf(stderr,"Parsed region for region %s\n", region);

    // need to seamlessly merge here with the dna2simplefeatures code
//    StringHash *readGroups = NULL;
//    readGroups = StringHash_new(STRINGHASH_SMALL);
//    for (i=0; i<Vector_getNumElement(intronBamConf->groupNames); i++) {
//      char *group = Vector_getElementAt(intronBamConf->groupNames, i);
//      if (verbosity > 0) fprintf(logfp, " %s", group);
//      StringHash_add(readGroups, group, &trueVal); 
//    }
//    if (verbosity > 0) fprintf(logfp, "\n");

    hts_itr_t *iter = sam_itr_queryi(RefineSolexaGenes_getAlignmentBamIndex(rsg), ref, begRange, endRange);

    fprintf(stderr,"before bam read loop... ");
    long counter[40];
    int j;
    for (j = 0; j < 40; j++) {
      counter[j] = 0;
    }
    while (bam_itr_next(rsg->alignment_sam, iter, read) >= 0) {
      if ((read->core.flag&BAM_FUNMAP) == 0) {
        uint8_t nm = bam_aux2A(bam_aux_get(read, "NM"));
        if (nm < 5) {
          if ((strand == 1 && !bam_is_rev(read)) || (strand == -1 && bam_is_rev(read))) {
  //          if (readGroups != NULL) { 
  //            uint8_t *rga = bam_aux_get(read,"RG");
  //            if (rga) {
  //              char *rg = bam_aux2Z(rga);
  //              if (!StringHash_contains(readGroups, rg)) {
  //                continue;
  //              }
  //            }
  //          }
            int32_t read_start = read->core.pos+1;
            int32_t read_end = bam_endpos(read)+1;
            int32_t real_start = 0;
            int32_t real_end = 40;
            if (read_start > pos) {
              real_start = read_start-pos+1;
            }
            if (real_end < pos+40) {
              real_end = read_end-pos+1;
            }
            for (int32_t i = real_start; i < real_end; i++) {
              ++counter[i];
            }
          }
        }
      }
    }
    fprintf(stderr,"Done\n");
    long max_count = counter[0];
    int max_index = 0;
    int low_index = 0;
    long low_count = counter[0];
    int i;
    fprintf(stderr, "%ld", counter[0]);
    for (i = 1; i < 60; i++) {
      fprintf(stderr, "\t%ld", counter[i]);
      if (counter[i] < low_count) {
        low_count = counter[i];
        low_index = i;
      }
      else if (counter[i] > max_count) {
        max_count = counter[i];
        max_index = i;
      }
    }
    fprintf(stderr, "\n");
    if (strand == 1) {
      lastpos = i+pos;
    }
    else {
      lastpos = pos-i;
    }
    if (low_index > max_index && low_index > 15) {
      if (low_count < (max_count*.5)) {
        printf(stderr, "Found end of 3'UTR at %ld low:%d max:%d\n", lastpos, low_count, max_count);
        if (strand == 1) {
          Exon_setEnd(exon, lastpos);
        }
        else {
          Exon_setStart(exon, lastpos);
        }
        bam_destroy1(read);
        sam_itr_destroy(iter);
        return lastpos;
      }
      else {
        printf(stderr, "Found possible end of 3'UTR at %ld low:%d max:%d\n", lastpos, low_count, max_count);
      }
    }
    else {
      printf(stderr, "Found possible end of 3'UTR at %ld low:%d max:%d\n", lastpos, low_count, max_count);
    }
    sam_itr_destroy(iter);
    res = strstr(res+6, "AATAAA");
  }
  bam_destroy1(read);
  return lastpos;
}


void RefineSolexaGenes_writeOutput(RefineSolexaGenes *rsg) {
  DBAdaptor *outdb = BaseGeneBuild_getDbAdaptor(rsg, RefineSolexaGenes_getOutputDb(rsg), 0, 0);

  int verbosity = RefineSolexaGenes_getVerbosity(rsg);

  GeneAdaptor *geneAdaptor = DBAdaptor_getGeneAdaptor(outdb);
// NIY:
//  $outdb->dbc->disconnect_when_inactive(0);

  Vector *output = RefineSolexaGenes_getOutput(rsg);
  
  int fails = 0;
  int total = 0;
//  GENE: 
  Analysis *anal = RefineSolexaGenes_getAnalysis(rsg);
  if (verbosity > 0) fprintf(stderr,"Analysis set to %s\n", Analysis_getLogicName(anal));

  int i;
  for (i=0; i<Vector_getNumElement(output); i++) {
    Gene *gene = Vector_getElementAt(output, i);

    Gene_setAnalysis(gene, anal);
    Gene_setSource(gene, RefineSolexaGenes_getBestScoreType(rsg));
    
    int j;
    for (j=0; j<Gene_getTranscriptCount(gene); j++) {
      Transcript *tran = Gene_getTranscriptAt(gene, j);
      Transcript_setAnalysis(tran, anal);
    }

    // filter single exon genes that may have been made through UTR trimming
    if (Gene_getExonCount(gene) == 1) {
      if (RefineSolexaGenes_getSingleExonModelType(rsg)) {
        Exon *exon = Transcript_getExonAt(Gene_getTranscriptAt(gene, 0), 0);
        if (Exon_getLength(exon) < RefineSolexaGenes_getMinSingleExonLength(rsg)) {
          if (verbosity > 0) fprintf(stderr, "Skipping writing single exon gene - length less than MinSingleExonLength\n");
          continue; // next GENE
        }
        Gene_setBiotype(gene, RefineSolexaGenes_getSingleExonModelType(rsg));
      } else {
        // dont store it
        if (verbosity > 0) fprintf(stderr, "Skipping writing single exon gene - no SingleExonModelType defined\n");
        continue;
      }
    }
    
// Was 'eval'ed
    GeneAdaptor_store(geneAdaptor, gene, 1);
/*
    if ($@){
      warning("Unable to store gene!!\n$@");
      $fails++;
    }
*/
    total++;
  }

// Can't be any currently, but maybe should be possible to trap errors
  if (fails > 0) {
    fprintf(stderr, "Not all genes could be written successfully (%d fails out of %d)\n", fails, total);
    exit(EX_IOERR);
  }

  if (verbosity > 0) fprintf(stderr, "Writing dna_align_features\n");
  if (RefineSolexaGenes_writeIntrons(rsg)) {
    Vector *intronFeatures = RefineSolexaGenes_getIntronFeatures(rsg);
  
    if (Vector_getNumElement(intronFeatures) > 1) {
      DNAAlignFeatureAdaptor *intronAdaptor = DBAdaptor_getDNAAlignFeatureAdaptor(outdb);
      fails = 0;
      total = 0;

      for (i=0; i<Vector_getNumElement(intronFeatures); i++) {
        DNAAlignFeature *intron = Vector_getElementAt(intronFeatures, i);

    // SMJS If want to edge match in apollo comment out these two lines
        DNAAlignFeature_setStart(intron, DNAAlignFeature_getStart(intron) + 1);
        DNAAlignFeature_setEnd(intron, DNAAlignFeature_getEnd(intron) - 1);
    //    eval {
    // tmpVec is just so can use _store method which takes a vector not a feature
    /* SMJS Perl did one at a time - I'm going to do all in one go after loop, so I don't have to make these temporary vectors
        Vector *tmpVec = Vector_new();
        Vector_addElement(tmpVec, intron);
        DNAAlignFeatureAdaptor_store(intronAdaptor, tmpVec);
        Vector_free(tmpVec);
    */
    //    };

    //    if ($@){
    //      warning("Unable to store DnaAlignFeature!!\n$@");
    //      $fails++;
    //    }
        total++;
      }

    // SMJS Moved from within loop - store all features in intronFeatures in one call
      DNAAlignFeatureAdaptor_store(intronAdaptor, intronFeatures);

      if (fails > 0) {
        fprintf(stderr, "Not all introns could be written successfully (%d fails out of %d)\n", fails, total);
        exit(EX_IOERR);
      }
    }
    else {
      fprintf(stderr, "Warning: No intron DNAAlignFeatures to write\n");
    }
  } else {
    fprintf(stderr, "NOTE: Not writing DNAAlignFeatures for introns\n");
  }
  if (total == 0) {
    fprintf(stderr, "No genes to write");
    exit(NO_GENES_WRITTEN);
  }
}


/*
=head2 ProcessTree

    Title        :   ProcessTree
    Usage        :   $self->ProcessTree
    Returns      :   String containing paths through the gene
    Args         :   A hash reference contianing the possible intron exons
                 :   Integer key for the hashref
                 :   String containing keys used up to this point
                 :   String containing the paths under construction
    Description  :   Recursive method that creates paths that explore all possible
                 :   routes through a hashref, uses a configurable recursion limit 
                 :   to prevent it running out of memory if the paths cannot be solved
                 :   or are too large to be practical

=cut
*/
/*
hashref contains lists of exon indexes keyed on intron hitseqname AND intron hitseqnames keyed on exon indexes. This allows alternating exon intron lookups with each ProcessTree
call
*/
/* 
Produces a hash of strings which have the format "exonIndex.intronName.exonIndex.intronName..."
*/

#define PT_ERROR 2
int RefineSolexaGenes_processTree(RefineSolexaGenes *rsg, StringHash *hashref, char *index, char *soFar, StringHash *paths, int *recLev) {
  int verbosity = RefineSolexaGenes_getVerbosity(rsg);
  (*recLev)++;
  // dont let it go on for ever eating memory
  if (limit > RefineSolexaGenes_getRecursiveLimit(rsg)){
    if (verbosity > 0) fprintf(stderr,"Too many recursive possibilities\n");
    (*recLev)--;
    return PT_ERROR;
  }

// Perl passes sofar by value so makes lots of strings with the different paths
// To replicate that, make copy of sofar here - must be allocated off heap, pointer to it is on stack, so should emulate the perl way
  sprintf(soFar,"%s%s.", soFar, index);
  //fprintf(stderr, "reclev = %d soFar is %s\n", *recLev, soFar);

  int nNode = 0;
  if (StringHash_contains(hashref, index)) {
    StringHash *subHash = StringHash_getValue(hashref, index);
    char **node = StringHash_getKeys(subHash);

    nNode = StringHash_getNumValues(subHash);

    char *savedSoFar;
    StrUtil_copyString(&savedSoFar, soFar, 0);

    int i;
    int failed = 0;
    for (i=0; i<nNode && !failed; i++) {
      char *child = node[i];

      limit++;
//  char tmpStr[512];
//  sprintf(tmpStr,"%s.",index);
//  newSoFar = StrUtil_appendString(newSoFar, tmpStr);
      //fprintf(stderr, "child = %s\n", child);
      int result = RefineSolexaGenes_processTree(rsg, hashref, child, soFar, paths, recLev);
      strcpy(soFar,savedSoFar);
      if (result == PT_ERROR) {
        limit = 0;
        (*recLev)--;
        failed = 1;
        //return PT_ERROR; 
      }
    }
    free(savedSoFar);
    for (i=0; i<nNode; i++) {
      free(node[i]);
    }
    free(node);
// Moved failed return to here so freeing of temporary data happens
    if (failed) return PT_ERROR;
  }
  if (nNode == 0) {
    //fprintf(stderr, "path key when nNode = 0 is %s\n", soFar);
    StringHash_add(paths, soFar, &trueVal);
  }
  (*recLev)--;
  return 0;
}

/*
=head2 ProcessTree

    Title        :   ProcessTree
    Usage        :   $self->process_paths
    Returns      :   String containing paths through the gene
    Args         :   A hash reference contianing the possible intron exons
                 :   Integer key for the hashref
                 :   String containing keys used up to this point
                 :   Integer flag indicating filtering should take place
    Description  :   Filters paths to remove the lowest scoring intron
                 :   for a given pair of exons where more than one intron
                 :   is possible. Filters progressivley if the paths cannot be
                 :   made for the model until the paths can be created or the 
                 :   model cannot be filtered any more, in this case the number
                 :   of recursions can be raised and the process repeated
                 :   untill the max_recursions limit is reached

=cut
*/
StringHash *RefineSolexaGenes_processPaths(RefineSolexaGenes *rsg, Vector *exons, Vector *exonIntron, StringHash *intronExon, int strict, int *giveUpFlag) {
  StringHash *variants = StringHash_new(STRINGHASH_SMALL);
  int removed = 0;
  int i;
  int j;
  int verbosity = RefineSolexaGenes_getVerbosity(rsg);

  //fprintf(stderr, "num element in exonIntron = %d\n", Vector_getNumElement(exonIntron));
  //fprintf(stderr, "num value in intronExon = %d\n", StringHash_getNumValues(intronExon));
  //# now lets make a hash of hashes holding which exons connect to which
  for (i=0; i<Vector_getNumElement(exons); i++) {
    Exon *exon = Vector_getElementAt(exons, i);
    if (i < Vector_getNumElement(exonIntron)) {
      Vector *exInti = Vector_getElementAt(exonIntron, i);
      char iStr[50];
      sprintf(iStr,"%d",i);
      if (exInti != NULL && Vector_getNumElement(exInti)) {
        if (strict) {
          //# Throw out exons that have retained introns for a start
  //        next if  $exons->[$i]->{'retained'};
          if (Exon_getFlags(exon) & RSGEXON_RETAINED) {
            continue;
          }
        }
  
        if (strict > 1) {
          // group all the introns by exon pairs
  // I think this code must have been moved from somewhere else. There doesn't seem to be a need for the intron_groups to be a hash because it only uses one key (i) to fill it (its local to the 'strict > 1' condition.
  // I've switched to using a single vector to store the intron group for this exon
          Vector *intronGroup = Vector_new();
          for (j=0; j<Vector_getNumElement(exInti); j++) {
            DNAAlignFeature *intron = Vector_getElementAt(exInti, j);
  
            if (strstr(DNAAlignFeature_getHitSeqName(intron), "REMOVED") == NULL) {
              Vector_addElement(intronGroup, intron);
            }
          }
  
          //# now lets sort these groups by score
          Vector_sort(intronGroup, SeqFeature_reverseScoreCompFunc);
  
          //# now lets see what they look like
          if (verbosity > 1) {
            fprintf(stderr,"EXON %d:%ld - %ld - %d\n",
                    i, Exon_getStart(exon), Exon_getEnd(exon), Exon_getStrand(exon));
          }
  
          for (j=0; j<Vector_getNumElement(intronGroup); j++) {
            DNAAlignFeature *intron = Vector_getElementAt(intronGroup, j);
            fprintf(stderr, "INTRON %d %s %ld %ld %d %.0f\n", i, DNAAlignFeature_getHitSeqName(intron), DNAAlignFeature_getStart(intron), DNAAlignFeature_getEnd(intron), DNAAlignFeature_getStrand(intron), DNAAlignFeature_getScore(intron));
          }
  
          if (Vector_getNumElement(intronGroup) > 1) {
            // remove the lowest scoring one
            //my $intron = pop( @{$intron_groups{$group}} ) ;
            DNAAlignFeature *intron = Vector_getLastElement(intronGroup);
            Vector_removeElementAt(intronGroup, Vector_getNumElement(intronGroup)-1);
  
            //fprintf(stderr,"Eliminating %s %.0f\n",DNAAlignFeature_getHitSeqName(intron), DNAAlignFeature_getScore(intron));
            char *hseqname = DNAAlignFeature_getHitSeqName(intron);
            if (!strstr(hseqname, "REMOVED")) {
  // Modifying the hseqname will mean it won't be found in intron_exon, so no paths will be generated using it
              char tmpStr[1024];
              sprintf(tmpStr, "%s-REMOVED", hseqname);
              DNAAlignFeature_setHitSeqName(intron, tmpStr);
              removed++;
            }
          }        
        }
        
        for (j=0; j<Vector_getNumElement(exInti); j++) {
          DNAAlignFeature *intron = Vector_getElementAt(exInti, j);
  
          if (StringHash_contains(intronExon, DNAAlignFeature_getHitSeqName(intron))) {
            //# only allow each intron to connect to 1 exon
            Vector *intEx = StringHash_getValue(intronExon, DNAAlignFeature_getHitSeqName(intron));
  
            int k;
            for (k=0; k<Vector_getNumElement(intEx); k++) {
              long exonInd = *((long *)Vector_getElementAt(intEx, k));
  
              char exonIndStr[1024];
              sprintf(exonIndStr,"%ld", exonInd);
  
              Exon *exon = Vector_getElementAt(exons, exonInd);
  // Note I inverted the condition
              if (DNAAlignFeature_getEnd(intron) <= Exon_getEnd(exon)) {
                //# store the possible paths as a hash (splice)variants
                char *hseqname = DNAAlignFeature_getHitSeqName(intron);
  
                fprintf(stderr,"adding variant for %s and %s %.0f\n", iStr, hseqname, DNAAlignFeature_getScore(intron));
                fprintf(stderr,"adding to variants %s hash value %s %.0f\n", iStr, hseqname, DNAAlignFeature_getScore(intron));
                //$variants->{$i}->{$intron->hseqname} = 1;
                if (!StringHash_contains(variants, iStr)) {
                  StringHash_add(variants, iStr, StringHash_new(STRINGHASH_SMALL));
                }
                fprintf(stderr,"adding to variants %s hash value %s %.0f\n",iStr, hseqname, DNAAlignFeature_getScore(intron));
                StringHash *varEHash = StringHash_getValue(variants, iStr);
                if (!StringHash_contains(varEHash, hseqname)) {
                  StringHash_add(varEHash, hseqname, &trueVal);
                } else {
                  fprintf(stderr,"Warning: varEHash already added variant %s to %s %.0f\n", hseqname, iStr, DNAAlignFeature_getScore(intron));
                }
  
                fprintf(stderr,"adding to variants %s hash value %s %.0f\n", hseqname, exonIndStr, DNAAlignFeature_getScore(intron));
                //$variants->{$intron->hseqname}->{$exon} = 1;
                if (!StringHash_contains(variants, hseqname)) {
                  StringHash_add(variants, hseqname, StringHash_new(STRINGHASH_SMALL));
                }
                StringHash *varIHash = StringHash_getValue(variants, hseqname);

                if (!StringHash_contains(varIHash, exonIndStr)) {
                  StringHash_add(varIHash, exonIndStr, &trueVal);
                } else {
                  fprintf(stderr,"Warning: varIHash already added variant %s to %s %.0f\n", exonIndStr, hseqname, DNAAlignFeature_getScore(intron));
                }
              }
            }
          }
        }
      }
    }
  }

  if (variants) {
    //fprintf(stderr, "Have %d values in variants\n", StringHash_getNumValues(variants));
  } else {
    if (verbosity > 0) fprintf(stderr, "variants is NULL\n");
  }
  
  if (strict && !removed) {
    Exon *firstExon = Vector_getElementAt(exons, 0);

    if (verbosity > 0) fprintf(stderr,"Warning: Cannot simplify this gene any more EXON 0: %ld - %ld - %d\n",
            Exon_getStart(firstExon), Exon_getEnd(firstExon), Exon_getStrand(firstExon));

    if (RefineSolexaGenes_getRecursiveLimit(rsg) < RefineSolexaGenes_getMaxRecursions(rsg)) {
      RefineSolexaGenes_setRecursiveLimit(rsg, RefineSolexaGenes_getRecursiveLimit(rsg) * 10);
      if (verbosity > 0) fprintf(stderr, "Upping recursive limit to %d to see if it helps\n", RefineSolexaGenes_getRecursiveLimit(rsg));
    } else {
      if (verbosity > 0) fprintf(stderr,"Giving up on EXON 0: %ld - %ld - %d\n",
              Exon_getStart(firstExon), Exon_getEnd(firstExon), Exon_getStrand(firstExon));
      StringHash_free(variants, StringHash_freeNoValFree);
//!!!!!!!!!!! NIY What to return
      *giveUpFlag = 1;
      return NULL;
    }
  }
  
  // work out all the possible paths given the features we have
  int result = 0;
  StringHash *paths = StringHash_new(STRINGHASH_LARGE);
  char *soFar = (char *)calloc(1000000,sizeof(char));
  for (i=0; i<Vector_getNumElement(exons) && result != PT_ERROR; i++) {
    limit = 0;
// Instead of undef for sofar I pass in ""
    char tmpStr[50];
    sprintf(tmpStr,"%d",i);
    int recLev = 0;
    soFar[0] = '\0';
    result = RefineSolexaGenes_processTree(rsg, variants, tmpStr, soFar, paths, &recLev);

    if (result == PT_ERROR) {
      if (verbosity > 0) fprintf(stderr, "Could not process cluster %d trying again with simpler cluster\n", i);
      //fprintf(stderr,"Freeing paths %p\n", paths);
      StringHash_free(paths, NULL);
      StringHash_free(variants, StringHash_freeNoValFree);
      return NULL;
    }
  }

  StringHash_free(variants, StringHash_freeNoValFree);

  free(soFar);
  return paths;
}

/*
=head2 model_cluster

    Title        :   model_cluster
    Usage        :   $self->model_cluster($models,$strand);
    Returns      :   2D array ref of exons and intron features
    Args         :   Array ref of  exons and intron features
                 :   Integer indicating strand
    Description  :   Clusters the initial models by start end
                 :   orders the models in each cluster by score

=cut
*/
Vector *RefineSolexaGenes_makeModelClusters(RefineSolexaGenes *rsg, Vector *models, int strand) {
  Vector *clusters = Vector_new();

  int verbosity = RefineSolexaGenes_getVerbosity(rsg);

  // sort them by the start of the fist exon ( the first array element )
  Vector_sort(models, Model_firstExonStartCompFunc);

  //my @models = sort { $a->[0]->start <=> $b->[0]->start }  @$models ;

  //# $model->[0]  = 1st exon 
  //# $model->[-3] = last exon 
  //# $model->[-2] = exon iuse
  //# $model->[-1] = score
  
  if (verbosity > 0) fprintf(logfp, "Have %d to cluster\n", Vector_getNumElement(models));

  int startInd = 0;

  int i;
  for (i=0; i<Vector_getNumElement(models); i++) {
    Model *model = Vector_getElementAt(models, i);
    int clustered = 0;

    int nClust = Vector_getNumElement(clusters);
    Exon *firstExon = Vector_getElementAt(model->features, 0);
    Exon *lastExon  = Vector_getLastElement(model->features);

    //fprintf(stderr, "Model start %ld end %ld\n",Exon_getStart(firstExon), Exon_getEnd(lastExon));
    int j;
    for (j=startInd; j<nClust; j++) {
//#    foreach my $cluster ( @clusters ) {
      // do they overlap?
      ModelCluster *cluster = Vector_getElementAt(clusters, j);

      if ( Exon_getStart(firstExon) <= cluster->end &&  Exon_getEnd(lastExon) >= cluster->start) {
        // Expand the cluster
        if (Exon_getStart(firstExon) < cluster->start) cluster->start = Exon_getStart(firstExon);
        if (Exon_getEnd(lastExon)    > cluster->end)   cluster->end   = Exon_getEnd(lastExon);

        Vector_addElement(cluster->models, model);
        clustered = 1;
      } else if (Exon_getStart(firstExon) > cluster->start) {
        startInd++;
      }
    }

    if (!clustered) {
      ModelCluster *cluster = ModelCluster_new();

// NIY: Do we need to add
      cluster->models = Vector_new();

      Vector_addElement(cluster->models, model);

      cluster->start  = Exon_getStart((Exon *)Vector_getElementAt(model->features,0));
      cluster->end    = Exon_getEnd((Exon *)Vector_getLastElement(model->features));
      cluster->strand = strand;

      Vector_addElement(clusters, cluster);
    }
  }

  if (verbosity > 0) fprintf(logfp, "Have %d after clustering them\n", Vector_getNumElement(clusters));
  return clusters;
}


/*
=head2 merge_exons

    Title        :   merge_exons
    Usage        :   $self->merge_exons($gene)
    Returns      :   Array ref of Bio::EnsEMBL::Exon
    Args         :   Bio::EnsEMBL::Gene
    Description  :   Merges adjacent exons where the intron is covered by repeats or
                 :   is very small

=cut
*/
// lets us merge exons with tiny  introns between them  unless they contain an intron
Vector *RefineSolexaGenes_mergeExons(RefineSolexaGenes *rsg, Gene *gene, int strand) {
  Vector *exons = Vector_new();

  if (Gene_getTranscriptCount(gene) < 1) {
// Perl was 'next' but that doesn't really make sense, should either be return or exit. I'll go with a warn and then return
    fprintf(stderr, "Warning: No transcript in gene\n");
    return exons;
  }

  Transcript *geneTrans = Gene_getTranscriptAt(gene, 0); 

  int i;
  for (i=0; i<Transcript_getExonCount(geneTrans); i++) {
    Exon *exon = Transcript_getExonAt(geneTrans, i);
    Vector_addElement(exons, ExonUtils_cloneExon(exon));
  }

  int ec = Vector_getNumElement(exons);

  StringHash *extraExons = RefineSolexaGenes_getExtraExons(rsg);
// For testing!!  StringHash *extraExons = StringHash_new(STRINGHASH_SMALL);

  // the extra exon is a list of start end coords of the spliced intron sections
  // ie: end:start:end:start where the 1st and last coords are anchors to tie it 
  // into our rough model both must match before we can try and add any potentialy 
  // novel exons in
  Vector *sortedExons = Vector_copy(exons);
  Vector_sort(sortedExons, SeqFeature_startCompFunc);

  Exon *lastExon = Vector_getLastElement(sortedExons);
  long endLastExon = Exon_getEnd(lastExon);
  Exon *firstExon = Vector_getElementAt(sortedExons, 0);
  long startFirstExon = Exon_getStart(firstExon);

  //char **keys = StringHash_getKeys(extraExons);
  char **keys = RefineSolexaGenes_getExtraExonsKeys(rsg);
  ExtraExonData **eeValues = RefineSolexaGenes_getExtraExonsValues(rsg);
  
  for (i=0; i<StringHash_getNumValues(extraExons); i++) {
    //ExtraExonData *eed = StringHash_getValue(extraExons, keys[i]);
    ExtraExonData *eed = eeValues[i];

//    my @coords = split(/:/,$key);

//    long startAnchor = shift(@coords);
//    long endAnchor = pop(@coords);
// NOTE: I'm not shifting and popping (ie removing from array) so need to take account of this when making exons

    long endAnchor   = eed->coords[eed->nCoord-1];
    if (startFirstExon > endAnchor) {
      continue;
    }

    long startAnchor = eed->coords[0];
    //#print "START AND END $start_anchor  $end_anchor \n";

    if (startAnchor > endAnchor) {
      fprintf(stderr, "START AND END %ld %ld\n", startAnchor,  endAnchor);
      exit(EX_SOFTWARE);
    }

    
    if (startAnchor > endLastExon) {
      break;
    }

    //# do the anchors lie within the model?
    //# SMJS Did some optimisation here
//#   foreach my $exon ( @exons ) {
//#      if ( $start_anchor <= $exon->end &&
//#           $start_anchor >= $exon->start ) {
//#        $start_anchor = -1;
//#      }
//#      if ( $end_anchor <= $exon->end &&
//#           $end_anchor >= $exon->start ) {
//#        $end_anchor = -1;
//#      }
//#    }

    if (RefineSolexaGenes_binSearchForOverlap(rsg, sortedExons, startAnchor)) {
      startAnchor = -1;
    }
    if (startAnchor == -1 && RefineSolexaGenes_binSearchForOverlap(rsg, sortedExons, endAnchor)) {
      endAnchor = -1;
    }
       
//#    my $i;
//#    for ($i=0;$i<scalar(@sortedexons); $i++) {
//#      my $exon = $sortedexons[$i];
//#      if ( $start_anchor <= $exon->end && 
//#           $start_anchor >= $exon->start ) {
//#        $start_anchor = -1;
//#        last;
//#      } elsif ($exon->start > $start_anchor) {
//#        last;
//#      }
//#    }
//#
//#    for ( ;$i<scalar(@sortedexons); $i++) {
//#      my $exon = $sortedexons[$i];
//#      if ( $end_anchor <= $exon->end && 
//#           $end_anchor >= $exon->start ) {
//#        $end_anchor = -1;
//#        last;
//#      } elsif ($exon->start > $end_anchor) {
//#        last;
//#      }
//#    }
    if ( startAnchor == -1 && endAnchor == -1 ) {
      // now to make new the exon(s)
// Note: Here start at loop 1 and end at nCoord-1 to avoid anchors (I haven't removed the anchors from array unlike perl)
      int j;
      for (j=1 ; j < eed->nCoord-1; j+=2) {
        long left  = eed->coords[j];
        long right = eed->coords[j+1];

        Exon *extra = RefineSolexaGenes_makeExon(rsg, left, right, eed->score, keys[i]);

        //my $extra = $self->make_exon(undef,$left,$right,$extra_exons->{$key},$key );

// Hack hack hack
        Exon_addFlag(extra, RSGEXON_EXTRA);

        Vector_addElement(exons, extra);
      }
    }
  }

/*
  for (i=0; i<StringHash_getNumValues(extraExons); i++) {
    free(keys[i]);
  }
  free(keys);
*/

// NIY: Free sortedExons
  Vector_free(sortedExons);

//  print "After extras add have " . scalar(@exons) . "\n";
  
//  print "Merging exons - done extras\n";

  
//#  @exons =  sort { $a->start <=> $b->start } @exons;
//#  # want to get rid of any overlapping exons
//#  while  ( $ec != scalar(@exons) ) {
//#    $ec = scalar(@exons);
//#    for ( my $i = 1 ; $i < scalar(@exons) ; $i++ ) {
//#      my $left_exon = $exons[$i-1];
//#      my $right_exon = $exons[$i];
//#      # do they overlap 
//#      if ( $left_exon->start <= $right_exon->end && 
//#           $left_exon->end >= $right_exon->start ) {
//#        # merge them 
//#        if (   $right_exon->end >= $left_exon->end &&
//#               $right_exon->start <= $left_exon->start ){
//#          print "HERE\n";
//#          $left_exon->{"_extra"} = 0;
//#        }
//#        $left_exon->start($right_exon->start) 
//#          if $right_exon->start < $left_exon->start;
//#        $left_exon->end($right_exon->end) 
//#          if $right_exon->end > $left_exon->end;
//#        # get rid of right exon
//#        splice(@exons,$i,1);
//#        $i-- ;
//#        @exons =  sort { $a->start <=> $b->start } @exons;
//#      }
//#    }
//#  }

  Vector_sort(exons, SeqFeature_startCompFunc);

// Note deliberately 1
  for (i=1 ; i<Vector_getNumElement(exons); i++) {
    Exon *leftExon  = Vector_getElementAt(exons, i-1);
    Exon *rightExon = Vector_getElementAt(exons, i);

     // do they overlap 
    if (Exon_getStart(leftExon) <= Exon_getEnd(rightExon) && 
        Exon_getEnd(leftExon)   >= Exon_getStart(rightExon)) {
      // merge them 
      if (Exon_getEnd(rightExon) >= Exon_getEnd(leftExon) &&
          Exon_getStart(rightExon) <= Exon_getStart(leftExon)) {
// Hack hack hack
        Exon_removeFlag(leftExon, RSGEXON_EXTRA);
      }
//# Don't see how this is possible - they are sorted on start!
//#      $left_exon->start($right_exon->start) 
//#        if $right_exon->start < $left_exon->start;
      
      if (Exon_getStart(rightExon) < Exon_getStart(leftExon)) {
        fprintf(stderr,"HOW IS THIS POSSIBLE\n");
        exit(EX_SOFTWARE);
      }

      if (Exon_getEnd(rightExon) > Exon_getEnd(leftExon)) {
        Exon_setEnd(leftExon, Exon_getEnd(rightExon)); 
      }

      // get rid of right exon
      //splice(@exons,$i,1);
      Exon *toRemove = Vector_removeElementAt(exons, i);
// NIY: Free removed exon??
      Exon_free(toRemove);

      i--;
//#      @exons =  sort { $a->start <=> $b->start } @exons;
    }
  }
  
//  fprintf(logfp, "Merging exons - done overlap filter\n");
  
  long offset = 0;
  for (i = 1; i<Vector_getNumElement(exons); i++) {
    Exon *exon     = Vector_getElementAt(exons, i);
    Exon *prevExon = Vector_getElementAt(exons, i-1);

    int intronCount = 0;

    Vector *introns = RefineSolexaGenes_fetchIntronFeatures(rsg, Exon_getEnd(prevExon), Exon_getStart(exon), &offset);
    
    // is the intron tiny?    
    // we know it lies across the boundary does it lie within the 2 exons?
    int j;
    for (j=0; j<Vector_getNumElement(introns); j++) {
      DNAAlignFeature *intron = Vector_getElementAt(introns, j);
 //#     print "INTRON " . $intron->start . " " . $intron->end . " " , $intron->strand ." " , $intron->score ."\n";

      //# ignore non consensus introns at this point
      if (DNAAlignFeature_getFlags(intron) & RSGINTRON_NONCANON) continue;
      //if (strstr(DNAAlignFeature_getHitSeqName(intron), "non canonical")) continue;

      if (DNAAlignFeature_getStart(intron) > Exon_getStart(prevExon) &&
          DNAAlignFeature_getEnd(intron) <  Exon_getEnd(exon) && 
          DNAAlignFeature_getStrand(intron) == strand) {
        intronCount++;
      }
    }

    // remove very small introns if there is no direct evidence for them
    if (Exon_getStart(exon) - Exon_getEnd(prevExon) <= 20  && intronCount == 0)   {
      Exon_setStart(exon, Exon_getStart(prevExon));

      //splice(@exons,$i-1,1);
      Exon *toRemove = Vector_removeElementAt(exons, i-1);
      Exon_free(toRemove);

      i--;
// Is this necessary???      next;
    }
    Vector_free(introns);
  }
//  fprintf(logfp, "Merging exons done\n");
  
  return exons;
}

Exon *RefineSolexaGenes_binSearchForOverlap(RefineSolexaGenes *rsg, Vector *exons, int pos) {
  int iMin = 0;
  int iMax = Vector_getNumElement(exons) - 1;

   while (iMax >= iMin) {
    int iMid = (iMax+iMin) / 2;

    Exon *e = Vector_getElementAt(exons, iMid);

    if (pos > Exon_getEnd(e)) {
      iMin = iMid + 1;
    } else if (pos < Exon_getStart(e)) {
      iMax = iMid - 1;
    } else {
      return e;
    }
  }

  return NULL;
}

/*
=head2 bam_2_intron_features
    Title        :   bam_2_intron_features
    Usage        :   $self->bam_2_intron_features($segment)
    Returns      :   None
    Args         :   Bam file segement
    Description  :   Fetches all alignments from the bam file segment, collapses them down into a 
                 :   non redundant set and builds a Bio::EnsEMBL::DnaDnaAlignFeature to 
                 :   represent it, then stores it in $self->intron_features
                 :   analyses splice sites for consensus and non consensus splices as this data is 
                 :   not stored in the BAM.
                 :   Also checks for small exons defined by a single read splicing at least twice
                 :   stores any additional exons found this way in $self->extra_exons

=cut
*/
void RefineSolexaGenes_bamToIntronFeatures(RefineSolexaGenes *rsg, IntronBamConfig *intronBamConf, htsFile *sam, bam_hdr_t *header, hts_idx_t *idx, int ref, int begRange, int endRange) {
  int verbosity = RefineSolexaGenes_getVerbosity(rsg);
  SliceAdaptor *sliceAdaptor = RefineSolexaGenes_getGeneSliceAdaptor(rsg);
  Vector *ifs = Vector_new();
  StringHash *extraExons = RefineSolexaGenes_getExtraExons(rsg);
  StringHash *idList = StringHash_new(STRINGHASH_LARGE);
  StringHash *readGroups = NULL;

  CigarBlock blockArray[1024];
  CigarBlock *mates[1024];
  int nMate;
  int i;
  for (i=0;i<1024;i++) {
    mates[i] = &blockArray[i];
  }


  if (intronBamConf->groupNames != NULL && Vector_getNumElement(intronBamConf->groupNames) > 0) {
    
    if (verbosity > 0) fprintf(logfp, "Limiting to read groups ");
    readGroups = StringHash_new(STRINGHASH_SMALL);
    for (i=0; i<Vector_getNumElement(intronBamConf->groupNames); i++) {
      char *group = Vector_getElementAt(intronBamConf->groupNames, i);
      if (verbosity > 0) fprintf(logfp, " %s", group);
      StringHash_add(readGroups, group, &trueVal); 
    }
    if (verbosity > 0) fprintf(logfp, "\n");
  }

  hts_itr_t *iter = sam_itr_queryi(idx, ref, begRange, endRange);
  bam1_t *read = bam_init1();

  int firstRead = 1;
  char name[1024];
  if (verbosity > 0) fprintf(stderr,"before bam read loop\n");
  while (bam_itr_next(sam, iter, read) >= 0) {

//READ:  
    char spliced;
    // ignore unspliced reads if the bam file is a mixture of spliced and 
    // unspliced reads
    if (intronBamConf->mixedBam) {
      uint8_t *spliceda = bam_aux_get(read, "XS");
      if (!spliceda) continue;
      spliced = bam_aux2A(spliceda);
    }
    // filter by read group if needed
    
    // need to recreate the ungapped features code as the
    // auto splitting code does not seem to work with > 2 features
    if (readGroups != NULL) { 
      uint8_t *rga = bam_aux_get(read,"RG");
      if (rga) {
        char *rg = bam_aux2Z(rga);
        if (!StringHash_contains(readGroups, rg)) {
          continue;
        }
      }
    }

    // Vector *mates = RefineSolexaGenes_getUngappedFeatures(rsg, read);
    nMate = RefineSolexaGenes_getUngappedFeatures(rsg, header, read, mates);

    //qsort(mates, nMate, sizeof(CigarBlock *), CigarBlock_startCompFunc);

    //Vector_sort(mates, CigarBlock_startCompFunc);
    //my @mates = sort { $a->[2] <=> $b->[2] } @{$self->ungapped_features($read)};

    // if mates > 2 then we have a possibility of adding in some extra exons into our rough models
    // as the read has spliced into and out of an exon
    // lets make them unique
    if (nMate > 2) {
      char keyString[2048]; keyString[0]='\0';
      long coords[1024]; // Hopefully won't have more than this!
      int nCoord = 0;

      for (i=0; i<nMate; i++) {
        CigarBlock *mate = mates[i];

        long start = mate->start;
        long end   = mate->end;

//  Unused      my $hstrand = $read->strand;
        if (i > 0) {
          sprintf(keyString, "%s%ld:", keyString, start);
          coords[nCoord++] = start;
        }
        if (i < nMate-1) {
          sprintf(keyString,"%s%ld:", keyString, end);
          coords[nCoord++] = end;
        }
      }

      ExtraExonData *eed = StringHash_getValue(extraExons, keyString);
      if (!eed) { //(!StringHash_contains(extraExons, keyString)) {
        eed = ExtraExonData_new(coords, nCoord);
      //  fprintf(stderr,"adding eed %s\n",keyString);
        StringHash_add(extraExons, keyString, eed);
      }
      //ExtraExonData *eed = StringHash_getValue(extraExons, keyString);

      eed->score++;
      //# print "Not doing extra_exon stuff for now\n";
    }

    //my $strand = $read->target->strand;
//    int strand = bam1_strand(read) == 0 ? -1 : 1;
    int strand = bam_is_rev(read) == 1 ? -1 : 1;
    if (intronBamConf->mixedBam) {
      if (spliced == '+') strand = 1;
      if (spliced == '-') strand = -1;
    } 

// Moved name setting out of loop
    if (firstRead == 1) {
      strcpy(name, header->target_name[read->core.tid]);
      StrUtil_strReplChr(name, '.', '*');
      //fprintf(stderr,"got target name %s\n", name);
      firstRead = 0;
    }

// Not used    long offset;
    char uniqueId[2048];
//    for (i=0; i<Vector_getNumElement(mates)-1; i++) {
//      CigarBlock *mate   = Vector_getElementAt(mates, i);
//      CigarBlock *mateP1 = Vector_getElementAt(mates, i+1);
    for (i=0; i<nMate-1; i++) {
      CigarBlock *mate   = mates[i];
      CigarBlock *mateP1 = mates[i+1];

      // intron reads should be split according to the CIGAR line
      // the default split function seems to ad
      // we want the ungapped features to make our introns
      // we dont allow . in the seq region name as we use them to delimit our paths
// Removed some unused vars here

      sprintf(uniqueId, "%s:%ld:%ld:%d", name, mate->end, mateP1->start, strand);

      IntronCoords *ic = StringHash_getValue(idList, uniqueId);
      if (!ic) {
        ic = IntronCoords_new(mate->end, mateP1->start, strand, -1, 0);
        StringHash_add(idList,
                       uniqueId,
                       ic);
      }
      //IntronCoords *ic = StringHash_getValue(idList, uniqueId);
      ic->score++;
    }
//    Vector_setFreeFunc(mates, CigarBlock_free);
//    Vector_free(mates);
  }
  sam_itr_destroy(iter);

  if (readGroups != NULL) StringHash_free(readGroups, NULL);

/* SMJS These are for testing
  // For param testing store introns with different anal to results
  my $conslim = $ENV{CONSLIM};
  my $nonconslim = $ENV{NONCONSLIM};
  my $intron_anal = $self->create_analysis_object("intron_c" . $conslim . "_nc" . $nonconslim);
*/

  if (verbosity > 0) fprintf(stderr,"before intron loop\n");
  //# collapse them down and make them into simple features
  IntronCoords **icArray = StringHash_getValues(idList);

  // May help smooth access to sequence to sort these - shouldn't be too costly to do the sort
  // May remove need for ifs sort below, but for now leave that
  qsort(icArray, StringHash_getNumValues(idList), sizeof(IntronCoords *), IntronCoords_intronStartCompFunc);

  Slice *chrSlice = RefineSolexaGenes_getChrSlice(rsg);
  char sliceName[2048];
  StrUtil_strReplChr(strcpy(sliceName, Slice_getName(chrSlice)), '.', '*');

  char sliceRegName[2048];
  StrUtil_strReplChr(strcpy(sliceRegName, Slice_getSeqRegionName(chrSlice)), '.', '*');

  Analysis *analysis = RefineSolexaGenes_getAnalysis(rsg);
  CachingSequenceAdaptor *cachingSeqAdaptor = DBAdaptor_getCachingSequenceAdaptor(Slice_getAdaptor(chrSlice)->dba);

  int min_intron_size = RefineSolexaGenes_getMinIntronSize(rsg)-1;
  for (i=0; i<StringHash_getNumValues(idList); i++) {
    IntronCoords *ic = icArray[i];

    // filter on score if appropriate
    if (intronBamConf->depth) {
      if (intronBamConf->depth > ic->score) {
        //#print "Rejecting on score " . $id_list{$key} ."\n";
        continue;
      }
    }
    long length =  ic->nextExonStart - ic->prevExonEnd -1;

    char name[2048];
    //I may loose some introns here but as they are smaller than the minimum intron size
    if (length > min_intron_size) {
      //sprintf(name,"%s:%ld:%ld:%d:", sliceName,
      sprintf(name,"%s:%ld:%ld:%d:", sliceRegName,
                                     ic->prevExonEnd+1,
                                     ic->nextExonStart-1,
                                     ic->strand);

      DNAAlignFeature *intFeat = DNAAlignFeature_new();

      DNAAlignFeature_setStart      (intFeat, ic->prevExonEnd);
      DNAAlignFeature_setEnd        (intFeat, ic->nextExonStart);
      DNAAlignFeature_setStrand     (intFeat, ic->strand);
      DNAAlignFeature_setHitStart   (intFeat, 1);
      DNAAlignFeature_setHitEnd     (intFeat, length);
      DNAAlignFeature_setHitStrand  (intFeat, 1);
      DNAAlignFeature_setSlice      (intFeat, chrSlice);
      DNAAlignFeature_setAnalysis   (intFeat, analysis);
      DNAAlignFeature_setScore      (intFeat, ic->score);
// Moved down      DNAAlignFeature_setHitSeqName (intFeat, name);
      char cigStr[256];
      sprintf(cigStr, "%ldM", length);
      DNAAlignFeature_setCigarString(intFeat, cigStr);

      int canonical = 1;
      // figure out if its cannonical or not
  // Was if->start+1 and if->start+2 rather than seqregionstart
/*
      Slice *leftSplice = SliceAdaptor_fetchByRegion(sliceAdaptor,
                                                     "toplevel",
                                                     Slice_getSeqRegionName(DNAAlignFeature_getSlice(intFeat)),
                                                     DNAAlignFeature_getSeqRegionStart(intFeat)+1,
                                                     DNAAlignFeature_getSeqRegionStart(intFeat)+2,
                                                     DNAAlignFeature_getSeqRegionStrand(intFeat),
                                                     NULL, 
                                                     0);
      Slice *rightSplice = SliceAdaptor_fetchByRegion(sliceAdaptor,
                                                      "toplevel",
                                                      Slice_getSeqRegionName(DNAAlignFeature_getSlice(intFeat)),
                                                      DNAAlignFeature_getSeqRegionEnd(intFeat)-2,
                                                      DNAAlignFeature_getSeqRegionEnd(intFeat)-1,
                                                      DNAAlignFeature_getSeqRegionStrand(intFeat),
                                                     NULL, 
                                                     0);

      //fprintf(stderr, "leftSplice seq = %s rightSplice seq = %s\n", Slice_getSeq(leftSplice), Slice_getSeq(rightSplice));
      char *leftSpliceSeq = Slice_getSeq(leftSplice);
      char *rightSpliceSeq = Slice_getSeq(rightSplice);
*/
      char *leftSpliceSeq = CachingSequenceAdaptor_fetchBySliceStartEndStrand(cachingSeqAdaptor, 
                                                                              chrSlice,
                                                                              DNAAlignFeature_getSeqRegionStart(intFeat)+1,
                                                                              DNAAlignFeature_getSeqRegionStart(intFeat)+2,
                                                                              DNAAlignFeature_getSeqRegionStrand(intFeat));
      char *rightSpliceSeq = CachingSequenceAdaptor_fetchBySliceStartEndStrand(cachingSeqAdaptor, 
                                                                               chrSlice,
                                                                               DNAAlignFeature_getSeqRegionEnd(intFeat)-2,
                                                                               DNAAlignFeature_getSeqRegionEnd(intFeat)-1,
                                                                               DNAAlignFeature_getSeqRegionStrand(intFeat));
      //fprintf(stderr,"%s %s\n", leftSpliceSeq, rightSpliceSeq);
      if (!strncmp(leftSpliceSeq, "NN", 2) && !strncmp(rightSpliceSeq, "NN", 2)) {
        fprintf(stderr,"Warning: Cannot find dna sequence for %s this is used in detecting non cannonical splices\n", name);
      } else {
        //# is it cannonical
        if (DNAAlignFeature_getStrand(intFeat) == 1 ) {
          // is it GTAG?
          if (strncmp(leftSpliceSeq, "GT", 2) || strncmp(rightSpliceSeq, "AG", 2)) {
            canonical = 0;
          }
        } else {
          //# is it GTAG?
          if (strncmp(rightSpliceSeq, "GT", 2) || strncmp(leftSpliceSeq, "AG", 2)) {
            canonical = 0;
          }
        }
      }
      free(leftSpliceSeq);
      free(rightSpliceSeq);
      //Slice_free(leftSplice);
      //Slice_free(rightSplice);

      char hitName[2048];
      if (canonical) {
        sprintf(hitName, "%s%s", name, "canon"); 
      } else {
        sprintf(hitName, "%s%s", name, "non canon"); 
        DNAAlignFeature_addFlag(intFeat, RSGINTRON_NONCANON);
      }
      DNAAlignFeature_setHitSeqName(intFeat, hitName);
  
      Vector_addElement(ifs, intFeat);
// NIY: Free splice site sequence strings and slices??
    }
  }

  StringHash_free(idList, IntronCoords_free);
  free(icArray);
  
  if (verbosity > 0) fprintf(stderr,"before filter\n");

  // sort them
  Vector_sort(ifs, SeqFeature_startCompFunc);
  if (RefineSolexaGenes_getFilterOnOverlapThreshold(rsg)) {
    Vector *tmpArray = Vector_new();
    int threshold = RefineSolexaGenes_getFilterOnOverlapThreshold(rsg);
    int is_one_threshold = RefineSolexaGenes_getIsOneThreshold(rsg);;

    int arrayLength = Vector_getNumElement(ifs);
    if (arrayLength > 1) {
      int j;
      for (j=0; j<arrayLength-1; j++) {
        DNAAlignFeature *ifj = Vector_getElementAt(ifs, j);
        if (DNAAlignFeature_getScore(ifj) > threshold) {
          Vector_addElement(tmpArray, ifj);
        }
        else {
          int k = 0;
          int overlappedSupport = 0;
          int is_one = DNAAlignFeature_getScore(ifj) <= is_one_threshold;
          int add = 0;
          int count = 0;
          while (add == 0) {
            ++k;
            DNAAlignFeature *ifjpk = Vector_getElementAt(ifs, j+k);

            if (((j+k) == arrayLength-1) || (DNAAlignFeature_getEnd(ifj) < DNAAlignFeature_getStart(ifjpk))) {
              //test the other side
              add = 2;
            }
            else {
              if (DNAAlignFeature_getStrand(ifj) != DNAAlignFeature_getStrand(ifjpk)) {
                continue;
              }
              if (DNAAlignFeature_getScore(ifj) < DNAAlignFeature_getScore(ifjpk)*0.05) {
                add = 1;
              }
              overlappedSupport += DNAAlignFeature_getScore(ifjpk);
              ++count;
            }
            //If the intron has only one read and it is overlapping intron with 50 reads or more, we skip it
            if (is_one && DNAAlignFeature_getScore(ifjpk) > 49) {
              add = 1;
            }
          }
          if (add == 2) {
            k = 0;
            while ((j-k) > 0) {
              ++k;
              DNAAlignFeature *ifjpk = Vector_getElementAt(ifs, j-k);

              if (DNAAlignFeature_getStrand(ifj) != DNAAlignFeature_getStrand(ifjpk)) {
                continue;
              }
              if ((DNAAlignFeature_getEnd(ifj) >= DNAAlignFeature_getStart(ifjpk)) &&
                (DNAAlignFeature_getStart(ifj) <= DNAAlignFeature_getEnd(ifjpk))) {
                if (DNAAlignFeature_getScore(ifj) < DNAAlignFeature_getScore(ifjpk)*0.05) {
                  add = 1;
                  k = j;
                }
                else {
                  overlappedSupport += DNAAlignFeature_getScore(ifjpk);
                  ++count;
                }
              }
              //If the intron has only one read and it is overlapping intron with 50 reads or more, we skip it
              if (is_one && DNAAlignFeature_getScore(ifjpk) > 49) {
                add = 1;
                k = j;
              }
            }
            if (add == 2) {
              Vector_addElement(tmpArray, ifj);
            }
          }
        }
      }
// NIY: Free old ifs??
      Vector_free(ifs);
      ifs = tmpArray;
    }
  }
//#  print STDERR 'RES: ', scalar(@ifs), "\n";

/* Filtering code is Steve's experimental code */
//# SMJS Filter here 

//  double nonConsLim = 15.0;
//  double consLim    = 5.0;
/*
  if (!defined($conslim) || !defined($nonconslim)) {
    die "Env vars for CONSLIM and NONCONSLIM not set\n";
  }
*/

  if (RefineSolexaGenes_getConsLim(rsg) > -0.1 && RefineSolexaGenes_getNonConsLim(rsg) > -0.1) {
    if (verbosity > 0) fprintf(stderr, "Filter parameters:  Consensus splice coverage %f  Non consensus splice coverage %f\n",
            RefineSolexaGenes_getConsLim(rsg), RefineSolexaGenes_getNonConsLim(rsg));
  
  // NIY: Need to free stuff here!!!!!
    int nCanon = 0;
    int nNonCanon = 0;
    Vector *tmp = Vector_new();
    for (i=0; i<Vector_getNumElement(ifs); i++) {
      DNAAlignFeature *f = Vector_getElementAt(ifs, i);
  
      if (DNAAlignFeature_getFlags(f) & RSGINTRON_NONCANON) {
        if (DNAAlignFeature_getScore(f) >= RefineSolexaGenes_getNonConsLim(rsg) && DNAAlignFeature_getLength(f) < 50000) {
          Vector_addElement(tmp, f);
          nNonCanon++;
        } else {
          //fprintf(stderr, "Rejected non_canonical feature with score %f\n", DNAAlignFeature_getScore(f));
          DNAAlignFeature_free(f);
        }
      } else {
        if (DNAAlignFeature_getScore(f) >= RefineSolexaGenes_getConsLim(rsg) && DNAAlignFeature_getLength(f) < 200000) {
          Vector_addElement(tmp, f);
          nCanon++;
        } else {
          //fprintf(stderr, "Rejected CANONICAL feature with score %f\n", DNAAlignFeature_getScore(f));
          DNAAlignFeature_free(f);
        }
      }
    }
    Vector_free(ifs);
    ifs = tmp;
  
    if (verbosity > 0) fprintf(stderr,"Got %d canonical introns and %d non canonical introns after filtering\n", nCanon, nNonCanon);
  }

  RefineSolexaGenes_setIntronFeatures(rsg, ifs);

// Do I need this - we fetch it at the start of the function and have been adding things to it, so it shouldn't need to be set
//  RefineSolexaGenes_setExtraExons(rsg, extraExons);

  if (verbosity > 0) fprintf(stderr,"Got %d unique introns  ", Vector_getNumElement(ifs));
  if (verbosity > 0) fprintf(stderr," and %d potential novel exons from %s\n", StringHash_getNumValues(extraExons), intronBamConf->fileName);
  return;
}


// Changed algorithm from perl - should have same outcome, but not require all the temporary objects
/* Aim:
     To get the match blocks at either end of every intron
   Algorithm:
     for each block
       If its a match
         Create block
         If we've just had an intron (hadIntron flag is set)
           unset currentBlock - we don't want to add a block twice
           add newly created block
           Unset hadIntron flag
         else
           if current block set at this point free it
           set current block to the new block
         end
       else If its an intron
         set hadIntron flag
         if current block is set
           add current block
         unset current block (so it doesn't get freed)
       else deal with pos offsets for other types
       end

     end for each block
*/
// NOTE: I think this method assumed that there were no cases with multiple 'M' blocks between two introns eg 20M 100N 20M D 20M 1000N 20M
//       I've added a check to make sure that assumption is true.
int RefineSolexaGenes_getUngappedFeatures(RefineSolexaGenes *rsg, bam_hdr_t *header, bam1_t *b, CigarBlock **ugfs) {
  CigarBlock *lastMatchBlock;
  int         hadIntron = 0;
  int         nIntron = 0;
  CigarBlock *currentBlock = NULL;
  uint32_t *  cigar = bam_get_cigar(b);
  int         nBlock = 0;
  CigarBlock  tmpBlock;

  int cigInd;
  int refPos;
  int readPos;
  for (cigInd = readPos = 0, refPos = b->core.pos+1; cigInd < b->core.n_cigar; ++cigInd) {
    int lenCigBlock = cigar[cigInd]>>4;
    int op          = cigar[cigInd]&0xf;

// M, =, X
    if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
      CigarBlock *block = CigarBlock_fill(&tmpBlock, CB_MATCH, refPos, refPos+lenCigBlock-1);

      if (hadIntron) {
        hadIntron = 0;
        CigarBlock_copy(ugfs[nBlock++], block);
        //Vector_addElement(ugfs, block);
        currentBlock = NULL;
      } else {
        currentBlock = block;
      }

      refPos += lenCigBlock; readPos += lenCigBlock;
// D
    } else if (op == BAM_CDEL) {
      refPos += lenCigBlock;
// S
    } else if (op == BAM_CSOFT_CLIP) {
      readPos += lenCigBlock;
// H
    } else if (op == BAM_CHARD_CLIP) {
// I
    } else if (op == BAM_CINS) {
      readPos += lenCigBlock;
// N
    } else if (op == BAM_CREF_SKIP) {
      hadIntron = 1;
      nIntron++;
      if (currentBlock) {
        //Vector_addElement(ugfs, currentBlock);
        CigarBlock_copy(ugfs[nBlock++], currentBlock);
        currentBlock = NULL;
      }
      refPos += lenCigBlock;
    } else {
      fprintf(stderr,"Cigar block type %d not supported\n", op);
      exit(EX_SOFTWARE);
    }
  }
  if (hadIntron) {
    kstring_t *bamline;
    if (sam_format1(header, b, bamline)) {
      fprintf(stderr, "Could not parse cigar line for %s\n", bam_get_qname(b));
      exit(EX_SOFTWARE);
    };
    fprintf(stderr,"Error parsing cigar string - don't have two M regions surrounding an N region (intron).\nBam entry = %s\n", ks_str(bamline));
    exit(EX_SOFTWARE);
  } 
  if (nBlock != nIntron+1  && nIntron > 0) {
    fprintf(stderr,"Error parsing cigar string - don't have the expected number of blocks (%d) for %d introns (have %d)\n", 
            nIntron+1, nIntron, nBlock);
    exit(EX_SOFTWARE);
  }

  return nBlock;

/* How perl did this:
  my @tmp_ugfs;
  my $string = $read->cigar_str;
  my $start = $read->start;
  my $end = $read->end;

  //#rint "THINGS $start $end $string\n";
  my @pieces = ( $string =~ /(\d*[MDN])/g );

  for my $piece ( @pieces ) {
    my ($length) = ( $piece =~ /^(\d*)/ );
    if( $length eq "" ) { $length = 1 }
    if( $piece =~ /M$/ ) {
      //
      // MATCH
      //
      my ( $qstart, $qend);
      $qstart = $start;
      $qend = $start + $length - 1;
      $start = $qend + 1;
      
      my $ugf;
      $ugf->[0] = $read->query->name;
      $ugf->[1] = $read->seq_id;
      $ugf->[2] = $qstart;
      $ugf->[3] = $qend;
      $ugf->[4] = $length."M";
      push @tmp_ugfs, $ugf;
      //#print "UNGAPPED " .$ugf->[2] .
        //#" " . $ugf->[3] . " " . $ugf->[4] ."\n";
    } elsif( $piece =~ /N$/ ) {
      //
      // INSERT
      //
      $start += $length;
      push @tmp_ugfs,"intron";
    } elsif( $piece =~ /D$/ ) {
      //
      // DELETION
      //
      $start += $length;
      push @tmp_ugfs,"deletion";
    } else {
      throw( "Illegal cigar line $string!" );
    }
  }

  // only return the UGFS either side of splices
  my %used_pieces;
  foreach ( my $i = 0 ; $i < scalar(@pieces); $i++ )  {
    my $piece = $pieces[$i];
    if ( $piece =~ /\d*N/) {
      // it's a splice push the Matches either side of it
      for ( my $j = $i-1 ; $j >= 0 ; $j-- ) {
        if ( $tmp_ugfs[$j] && $pieces[$j] =~ /\d*M/ )  {
          unless ( $used_pieces{$j} ) {
            my $ugf =  $tmp_ugfs[$j];
            $self->throw("Cannot find ugf $j\n") unless $ugf;
            push @ugfs, $ugf;
            $used_pieces{$j} =1;
            last;
          }
        }
      }

      for ( my $j = $i+1 ; $j < scalar(@pieces)  ; $j++ ) {
        if ( $tmp_ugfs[$j] && $pieces[$j] =~ /\d*M/ )  {
          unless ( $used_pieces{$j} ) {
            my $ugf =  $tmp_ugfs[$j];
            $self->throw("Cannot find ugf $j\n") unless $ugf;
            push @ugfs, $ugf;
            $used_pieces{$j} =1;
            last;
          }
        }
      }
    }
  }

  return $ugfs;
*/
}


/*
=head2 dna_2_extra_exons
    Title        :   dna_2_extra_exons
    Usage        :   $self->dna_2_extra_exons($start,$end)
    Returns      :   None
    Args         :   Int start
                 :   Int end
    Description  :   Fetches all dna_align_features from the small intron db that lie within
                 :   the range determined by start and end, collapses them down into a 
                 :   non redundant set and builds extra exons to account for short exons
                 :   missed by genomic alignemnt of reads

=cut
*/

/* Was unused in perl (never called), and is not implemented correctly currently, because it makes an array of exons in extraExons whereas bam method makes a hash of start:end's

Vector *RefineSolexaGenes_dnaToExtraExons(RefineSolexaGenes *rsg, long start, long end) {

  Vector *extraIntrons = Vector_getIntronFeatures(rsg);

  SliceAdaptor *smallIntronSliceAdaptor = RefineSolexaGenes_getSmallIntronSliceAdaptor(rsg);
  Slice *smallIntronSlice = SliceAdaptor_fetchByRegion(smallIntronSliceAdaptor, "toplevel", Slice_getSeqRegionName(chrSlice), start, end); 

  my @ifs;
  push @ifs, @$extra_introns if $extra_introns;
  my $rough_genes = $self->prelim_genes;
  my %exon_list;
  my %id_list;
  
  // fetch all the dna_align_features for this slice by logic name
  my @extra_reads;
  // look for extra introns from ESLA

  // featch all the dna_align_features for this slice by logic name
  Vector *logicNames = RefineSolexaGenes_getLogicNames(rsg);
  Vector *reads;

  if (Vector_getNumElement(logicNames)) {
    reads = Vector_new();
    fprintf(stderr,"Fetching reads with logic names: ");

    int i;
    for (i; i<Vector_getNumElement(logicNames); i++) {
      char *logicName = Vector_getElementAt(logicNames, i);

      fprintf(stderr,"%s ", logicName);
      Vector *feats = Slice_getAllDNAAlignFeatures(smallIntronSlice, logicName);
      Vector_append(reads, feats); 
      Vector_free(feats);
    }
    fprintf(stderr,"\n");
  } else {
  //# fetch them all if no logic name is Supplied
    reads =  Slice_getAllDNAAlignFeatures(smallIntronSlice);
  }
  fprintf(stderr, "Got %d extra reads\n", Vector_getNumElement(reads));

  // process extra reads and assign them to rough models
  while ( scalar @extra_reads > 0 ) {
    char *roughId = NULL;

    my $read = pop(@extra_reads);
    my $type = 'canonical';
    $type = 'non canonical' if ( $read->hseqname =~ /\:NC$/ ) ;
    $read = $read->transfer($self->chr_slice);
    // assign a rough model to the reads
    my @ugfs = $read->ungapped_features;
    @ugfs = sort { $a->start <=> $b->start } @ugfs;
    next if ( scalar(@ugfs) == 0 ) ;
    // first and or last ugf should overlap an exon
  ROUGH:  foreach my $rough ( @$rough_genes ) {
      foreach my $exon ( @{$rough->get_all_Exons} ) {
        if( ( $exon->start <= $ugfs[0]->end &&  $exon->end >= $ugfs[0]->start ) or 
            ( $exon->start <= $ugfs[-1]->end &&  $exon->end >= $ugfs[-1]->start ) ) {
          //ONLY USE THIS READ WITH THIS MODEL
          $roughid = $rough->stable_id;
          last ROUGH;
        }
      }
    }
    unless ( $roughid ) {
      print STDERR "Ignoring read " . $read->hseqname . " cannot pair it with a rough model\n";
    }
    for ( my $i = 0 ; $i < scalar(@ugfs) - 1 ; $i++ ) {
      //# one read can span several exons so make all the features 
      // cache them by internal boundaries
      // we use . to deliminate entries in our paths so dont allow them in the seq_region_name or it wont work
      char name[1024];
      strcpy(name, DNAAlignFeature_getSeqRegionName(read));
      StrUtil_strReplChr(name, '.', '*');
      
      char uniqueId[2048];
      sprintf(uniqueId, "%s:%ld:%ld:%d:%s", name, DNAAlignFrature_getEnd(ugf), 
              DNAAlignFeature_getStart(ugfP1), DNAAlignFeature_getStrand(read), type)

      $id_list{$unique_id} ++;
      if  ( $i > 0 && $i < scalar(@ugfs) - 1 ) {
        // if the read splices into and out of an exon, store that exon in case it is not 
        // found in the rough model but is an extra short exon determined by ExonerateSolexaLocalAlignment
        // check that it is introns on both sides and not inserts, need to be bigger than the min intron length
        next unless $ugfs[$i+1]->start - $ugfs[$i]->end > $self->MIN_INTRON_SIZE;
        next unless $ugfs[$i]->start - $ugfs[$i-1]->end > $self->MIN_INTRON_SIZE;

        strcpy(name, DNAAlignFeature_getSeqRegionName(read));
        StrUtil_strReplChr(name, '.', '*');
        
        char uniqueId[2048];
        sprintf(uniqueId, "%s:%ld:%ld:%d:%s", name, DNAAlignFrature_getStart(ugf), 
                DNAAlignFeature_getEnd(ugf), DNAAlignFeature_getStrand(read), type)
        $exon_list{$roughid}{$unique_id} = $ugfs[$i];        
      }
    }
  }
  
  //collapse down the intron feaures and make them into simple features
  foreach my $key ( keys %id_list ) {
   // print "KEY $key\n";
    my @data = split(/:/,$key) ;
    my $length =  $data[2] - $data[1] -1;
    next unless $length > 0 ;
    my $name = $self->chr_slice->seq_region_name . ":" . ($data[1]+1) . ":" . ($data[2] -1) .":" . $data[3] .":".$data[4];
    my $if = Bio::EnsEMBL::DnaDnaAlignFeature->new
      (
       -start => $data[1],
       -end => $data[2],
       -strand => $data[3],
       -hstart => 1,
       -hend => $length,
       -hstrand => 1,
       -slice => $self->chr_slice,
       -analysis => $self->analysis,
       -score =>  $id_list{$key},
       -hseqname => $name,
       -cigar_string => $length ."M",
      );
    push @ifs , $if;
  }

  // sort them
  @ifs = sort {$a->start <=> $b->start} @ifs;
  $self->intron_features(\@ifs);
  print STDERR "Got " . scalar(@ifs) . " intron features\n";
  // make the exon features
  my @extra_exons;
  foreach my $rough ( keys %exon_list ) {
    foreach my $key ( keys %{$exon_list{$rough}} ) {
      my $extra_exon = $self->make_exon($exon_list{$rough}{$key});

// Hack hack hack
      $extra_exon->{"_extra"} = 1;
      $extra_exon->{"_model"} = $rough;

      Vector_addElement(push @extra_exons, $extra_exon;
    }
  }

  fprintf(stderr, "Got %d extra exons\n", Vector_getNumElement(extraExons));

  // make the exon features
  RefineSolexaGenes_setExtraExons(rsg, extraExons);
  return;
}
*/


/*
=head2 dna_2_intron_features
    Title        :   dna_2_intron_features
    Usage        :   $self->dna_2_intron_features($start,$end)
    Returns      :   None
    Args         :   Int start
                 :   Int end
    Description  :   Fetches all dna_align_features from the intron db that lie within
                 :   the range determined by start and end, collapses them down into a 
                 :   non redundant set and builds a Bio::EnsEMBL::DnaAlignFeature to 
                 :   represent it, then stores it in $self->intron_features
                 :   also checks for small exons defined by a single read splicing at least twice
                 :   stores any additional exons found this way in $self->extra_exons

=cut
*/
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// Note NOTE Note: This code is not currently in use and probably needs some fixing up to match the Bam reading code
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
void RefineSolexaGenes_dnaToIntronFeatures(RefineSolexaGenes *rsg, long start, long end) {
  int verbosity = RefineSolexaGenes_getVerbosity(rsg);
// Unused  my $rough_genes = $self->prelim_genes;
  SliceAdaptor *intronSliceAdaptor = RefineSolexaGenes_getIntronSliceAdaptor(rsg);
  Slice *chrSlice = RefineSolexaGenes_getChrSlice(rsg);

  Slice *intronSlice = SliceAdaptor_fetchByRegion(intronSliceAdaptor, "toplevel", Slice_getSeqRegionName(chrSlice), start, end, 1, NULL, 0); 

  // fetch all the dna_align_features for this slice by logic name
  Vector *logicNames = RefineSolexaGenes_getLogicNames(rsg);
  Vector *reads;

  if (logicNames != NULL && Vector_getNumElement(logicNames)) {
    reads = Vector_new();
    if (verbosity > 0) fprintf(stderr,"Fetching reads with logic names: ");

    int i;
    for (i=0; i<Vector_getNumElement(logicNames); i++) {
      char *logicName = Vector_getElementAt(logicNames, i);

      if (verbosity > 0) fprintf(stderr,"%s ", logicName);
      Vector *feats = Slice_getAllDNAAlignFeatures(intronSlice, logicName, NULL, NULL, NULL);
      Vector_append(reads, feats); 
      Vector_free(feats);
    }
    if (verbosity > 0) fprintf(stderr,"\n");
  } else {
  //# fetch them all if no logic name is Supplied
    reads =  Slice_getAllDNAAlignFeatures(intronSlice, NULL, NULL, NULL, NULL);
  }
  if (verbosity > 0) fprintf(stderr, "Got %d reads\n", Vector_getNumElement(reads));

  StringHash *idList = StringHash_new(STRINGHASH_MEDIUM);

  int i;
  for (i=0; i<Vector_getNumElement(reads); i++) {
    DNAAlignFeature *read = Vector_getElementAt(reads, i);

    int isCanonical = 1;
 // Was  $type = 'non canonical' if ( $read->hseqname =~ /\:NC$/ ) ;
    if (strstr(DNAAlignFeature_getHitSeqName(read), ":NC")) {
      isCanonical = 0;
    }
    
    // Do we actually need to tranfer this???
    DNAAlignFeature *origRead = read;
// Use SeqFeature_transfer as it should be all that's needed
    read = SeqFeature_transfer((SeqFeature *)read, chrSlice);
// NIY: Need to free

    Vector *ugfs = DNAAlignFeature_getUngappedFeatures(read);

    if (Vector_getNumElement(ugfs) > 0) {
      Vector_sort(ugfs, SeqFeature_startCompFunc);
     
      //# one read can span several exons so make all the features 
      // cache them by internal boundaries
      // we use . to deliminate entries in our paths so dont allow them in the seq_region_name or it wont work
      char name[1024];
      strcpy(name, Slice_getSeqRegionName(DNAAlignFeature_getSlice(read)));
      StrUtil_strReplChr(name, '.', '*');
      
      char uniqueId[2048];

      int j;
      for (j=0 ; j<Vector_getNumElement(ugfs)-1; j++) {
        DNAAlignFeature *ugf   = Vector_getElementAt(ugfs, j);
        DNAAlignFeature *ugfP1 = Vector_getElementAt(ugfs, j+1);
        
        sprintf(uniqueId, "%s:%ld:%ld:%d:%s", name, DNAAlignFeature_getEnd(ugf), 
                DNAAlignFeature_getStart(ugfP1), DNAAlignFeature_getStrand(read), isCanonical ? "canon" : "non canon");

        if (!StringHash_contains(idList, uniqueId)) {
          StringHash_add(idList, 
                         uniqueId, 
                         IntronCoords_new(DNAAlignFeature_getEnd(ugf), DNAAlignFeature_getStart(ugfP1), DNAAlignFeature_getStrand(read), isCanonical, 0));
        }
        IntronCoords *ic = StringHash_getValue(idList, uniqueId);
        ic->score++;
      }
    }
    Vector_setFreeFunc(ugfs,DNAAlignFeature_freeImpl);
    Vector_free(ugfs);
  }

  if (verbosity > 0) fprintf(stderr, "Got %d collapsed introns\n", StringHash_getNumValues(idList));

  //# collapse them down and make them into dna align features


  Vector *intFeats = Vector_new();
  Analysis *analysis = RefineSolexaGenes_getAnalysis(rsg);

  char sliceName[2048];
  StrUtil_strReplChr(strcpy(sliceName, Slice_getName(chrSlice)), '.', '*');

  char sliceRegName[2048];
  StrUtil_strReplChr(strcpy(sliceRegName, Slice_getSeqRegionName(chrSlice)), '.', '*');

  IntronCoords **icArray = StringHash_getValues(idList);

  for (i=0; i<StringHash_getNumValues(idList); i++) {
    IntronCoords *ic = icArray[i];
    long length =  ic->nextExonStart - ic->prevExonEnd -1;

    char name[2048];
    if (length > 0) {
      sprintf(name,"%s:%ld:%ld:%d:%s", sliceRegName, 
                                       ic->prevExonEnd+1, 
                                       ic->nextExonStart-1, 
                                       ic->strand, 
                                       ic->isCanonical == 0 ? "non canon" : "canon");

      DNAAlignFeature *intFeat = DNAAlignFeature_new();

      DNAAlignFeature_setStart      (intFeat, ic->prevExonEnd);
      DNAAlignFeature_setEnd        (intFeat, ic->nextExonStart);
      DNAAlignFeature_setStrand     (intFeat, ic->strand);
      DNAAlignFeature_setHitStart   (intFeat, 1);
      DNAAlignFeature_setHitEnd     (intFeat, length);
      DNAAlignFeature_setHitStrand  (intFeat, 1);
      DNAAlignFeature_setSlice      (intFeat, chrSlice);
      DNAAlignFeature_setAnalysis   (intFeat, analysis);
      DNAAlignFeature_setScore      (intFeat, ic->score);
      DNAAlignFeature_setHitSeqName (intFeat, name);
      char cigStr[256];
      sprintf(cigStr, "%ldM", length);
      DNAAlignFeature_setCigarString(intFeat, cigStr);

      if (ic->isCanonical == 0) {
        DNAAlignFeature_addFlag(intFeat, RSGINTRON_NONCANON);
      }

      Vector_addElement(intFeats, intFeat);
    }
  }

  free(icArray);

  StringHash_free(idList, IntronCoords_free);

  // sort them
  Vector_sort(intFeats, SeqFeature_startCompFunc);

  RefineSolexaGenes_setIntronFeatures(rsg, intFeats);

  if (verbosity > 0) fprintf(stderr, "Got %d intron features\n", Vector_getNumElement(intFeats));
  return;
}

/*
=head2 fetch_intron_features

    Title        :   fetch_intron_features
    Usage        :   $self->fetch_intron_features($start,$end)
    Returns      :   Array ref of Bio::EnsEMBL::DnaAlignFeature
    Args         :   Int start
                 :   Int end
    Description  :   Accesses the pre computed simple features representing introns
                 :   Filters out non consensus models that overlap consensus models

=cut
*/

// Not neat - not trying to find an exact boundary, just to find an approximate start point in the array which
// is before the wanted start point - the max intron length in the vector.
// aim is to reduce the number of comparisons done in fetchIntronFeatures loop
int RefineSolexaGenes_binSearchIntrons(RefineSolexaGenes *rsg, Vector *introns, long pos) {
  int iMin = 0;
  int iMax = Vector_getNumElement(introns) - 1;
  int nLoop = 0;

  long maxLength = RefineSolexaGenes_getLongestIntronLength(rsg);

  pos = pos - maxLength - 1; // -1 is just to make sure

  while (iMax >= iMin) {
    int iMid = (iMax+iMin) / 2;
    nLoop++;

    DNAAlignFeature *daf = Vector_getElementAt(introns, iMid);

    if (pos > DNAAlignFeature_getStart(daf)) {
      iMin = iMid + 1;
// Was <= but don't think it should be, so switched to <
    } else if (pos < DNAAlignFeature_getStart(daf)) {
      iMax = iMid - 1;
    } else {
      break;
    }
  }

  int indToReturn = iMin < iMax ? iMin : iMax;
  
  int nSecondLoop = 0;
  if (indToReturn > 0) {
    DNAAlignFeature *daf = Vector_getElementAt(introns, indToReturn);
    while (indToReturn > 0 && DNAAlignFeature_getStart(daf) > pos) {
      //fprintf(stderr," feature start = %ld pos = %ld\n", DNAAlignFeature_getStart(daf) , pos);
      daf = Vector_getElementAt(introns, --indToReturn);
      nSecondLoop++;
    }
  } else if (indToReturn < 0) {
    indToReturn = 0;
  }

  //fprintf(stderr,"nLoop in binsearch = %d searching for pos %ld nSecondLoop %d\n",nLoop, pos, nSecondLoop);
  //fprintf(stderr,"  ended with iMin of %d and iMax of %d\n", iMin, iMax);
  //DNAAlignFeature *daf = Vector_getElementAt(introns, indToReturn);
  //fprintf(stderr,"  feature at index %d start %ld end %ld\n", indToReturn, DNAAlignFeature_getStart(daf), DNAAlignFeature_getEnd(daf));
  
  return indToReturn;
}

Vector *RefineSolexaGenes_fetchIntronFeatures(RefineSolexaGenes *rsg, long start, long end, long *offsetP) {
  Vector *sfs = RefineSolexaGenes_getIntronFeatures(rsg);

  Vector *filteredIntrons = Vector_new();
  if (Vector_getNumElement(sfs)) {
    long intronStart = 0;
    if (*offsetP) {
      intronStart = *offsetP;
    } else {
      intronStart = RefineSolexaGenes_binSearchIntrons(rsg, sfs, start);
    }

    int index = -1;

    Vector *chosenSf = Vector_new();

    Vector_setBatchSize(chosenSf, 200);

    //fprintf(stderr, "Have %d intron features in fetchIntronFeatures. OffsetP = %ld intronStart = %ld\n", Vector_getNumElement(sfs), *offsetP, intronStart);

    // sfs is a sorted array
    int i;
    for (i=intronStart; i<Vector_getNumElement(sfs); i++) {
      DNAAlignFeature *intron = Vector_getElementAt(sfs, i);

  /*
      fprintf(stderr, " in first loop %ld %ld %d %f %s\n",
              DNAAlignFeature_getStart(intron),
              DNAAlignFeature_getEnd(intron),
              DNAAlignFeature_getStrand(intron),
              DNAAlignFeature_getScore(intron),
              DNAAlignFeature_getHitSeqName(intron));
  */

      if (DNAAlignFeature_getStart(intron) > end) break;

      if (DNAAlignFeature_getStart(intron) <= end && DNAAlignFeature_getEnd(intron) >= start) {
        Vector_addElement(chosenSf, intron);

        // remember the position of the 1st intron to overlap
        // this exon - we will start counting from here next time
  // Note perl used unless $index so 0 is ignored - I've put in index==0 to emulate that but its a bit horrid and probably a bug
  // in the perl
  //      if (index == -1 || index == 0) {
  // Do as <1 for now (hopefully faster) - but take note of above note about perl behaviour
        if (index < 1) {
   //       fprintf(stderr,"Setting index to %d\n", index);
          index = i;
        }
      }
    }

  // INTRON:
    int startCompInd = 0;
    for (i=0; i<Vector_getNumElement(chosenSf); i++) {
      DNAAlignFeature *intron = Vector_getElementAt(chosenSf, i);

  /*
      fprintf(stderr, " in loop %ld %ld %d %f %s\n",
              DNAAlignFeature_getStart(intron),
              DNAAlignFeature_getEnd(intron),
              DNAAlignFeature_getStrand(intron),
              DNAAlignFeature_getScore(intron),
              DNAAlignFeature_getHitSeqName(intron));
  */

  //    if (strstr(DNAAlignFeature_getHitSeqName(intron), "non canonical"))
      if (DNAAlignFeature_getFlags(intron) & RSGINTRON_NONCANON) {
        // check it has no overlap with any consensus introns
        // unless it out scores a consensus intron
        int done = 0;
        int j;
        for (j=startCompInd; j<Vector_getNumElement(chosenSf) && !done; j++) {
          DNAAlignFeature *compIntron = Vector_getElementAt(chosenSf, j);

          // Optimisation - break out of loop once there can no longer be overlap (start of compIntron > end of intron) - remember introns are sorted
          // Doesn't seem to help performance much in tests but I think maybe in extreme cases it will help so leave for now
          if (DNAAlignFeature_getEnd(intron) < DNAAlignFeature_getStart(compIntron)) break;

          //if (strstr(DNAAlignFeature_getHitSeqName(compIntron), "non canonical") == NULL)
          if (!(DNAAlignFeature_getFlags(compIntron) & RSGINTRON_NONCANON)) {
            if (DNAAlignFeature_getStrand(intron) == DNAAlignFeature_getStrand(compIntron) &&
                DNAAlignFeature_getEnd(intron)    >  DNAAlignFeature_getStart(compIntron) &&
                DNAAlignFeature_getStart(intron)  <  DNAAlignFeature_getEnd(compIntron)) {
              // Done flag replaces thisnext INTRON if $intron->score <= $i->score;
              if (DNAAlignFeature_getScore(intron) <= DNAAlignFeature_getScore(compIntron)) {
                done = 1;
              }
            }
          }
        }
        if (!done) {
          Vector_addElement(filteredIntrons, intron);
        }
      } else {
        Vector_addElement(filteredIntrons, intron);
      }
    }

    //fprintf(stderr,"index = %d\n",index);
    if (index != -1) {
      *offsetP = index;
    }


    //fprintf(stderr,"filteredIntrons :\n");
  /*
    for (i=0; i<Vector_getNumElement(filteredIntrons); i++) {
      DNAAlignFeature *intron = Vector_getElementAt(filteredIntrons, i);
      fprintf(stderr, " %ld %ld %d %f %s\n",
              DNAAlignFeature_getStart(intron),
              DNAAlignFeature_getEnd(intron),
              DNAAlignFeature_getStrand(intron),
              DNAAlignFeature_getScore(intron),
              DNAAlignFeature_getHitSeqName(intron));
    }
  */

    Vector_free(chosenSf);
  }
  
  return filteredIntrons;
}


/*

=head2 make_exon
    Title        :   pad_exons
    Usage        :   $self->($ungapped_feature)
    Returns      :   Bio::EnsEMBL::Exon
    Args         :   Bio::EnsEMBL::FeaturePair 
    Description  :   Takes an ungapped feature, pads it and builds a 
                 :   Exon from it 
=cut
*/
// Note : Removed ability to call with a ugf which was only used in unused method
Exon *RefineSolexaGenes_makeExon(RefineSolexaGenes *rsg, long start, long end, double score, char *displayId) {
  Slice *chrSlice = RefineSolexaGenes_getChrSlice(rsg);

/*
  if (ugf) {
    start     = ugf->start;
    end       = ugf->end;
    displayId = ugf->displayId;
    score     = ugf->score;
  }
*/

  Exon *paddedExon = ExonUtils_createExon(start-20, end+20, -1, -1, -1, RefineSolexaGenes_getAnalysis(rsg), NULL, 0, RefineSolexaGenes_getChrSlice(rsg), NULL, 0);


  // dont let it fall of the slice because of padding
  if (Exon_getStart(paddedExon) <= 0) {
    Exon_setStart(paddedExon, 1);
  }

  if (Exon_getEnd(paddedExon) >= Slice_getLength(chrSlice)) {
    Exon_setEnd(paddedExon, Slice_getLength(chrSlice)-1);
  }
  
  DNAAlignFeature *feat = DNAAlignFeature_new();
  DNAAlignFeature_setSlice     (feat, RefineSolexaGenes_getChrSlice(rsg));
  DNAAlignFeature_setStart     (feat, Exon_getStart(paddedExon));
  DNAAlignFeature_setEnd       (feat, Exon_getEnd(paddedExon));
  DNAAlignFeature_setStrand    (feat, -1); // Why -1??????????????
  DNAAlignFeature_setHitSeqName(feat, displayId);
  DNAAlignFeature_setHitStart  (feat, 1);
  DNAAlignFeature_setHitStrand (feat, 1);
  DNAAlignFeature_setHitEnd    (feat, Exon_getLength(paddedExon));
  DNAAlignFeature_setAnalysis  (feat, RefineSolexaGenes_getAnalysis(rsg));
  DNAAlignFeature_setScore     (feat, score);

  char tmpStr[256];
  sprintf(tmpStr, "%ldM", Exon_getLength(paddedExon));

// Only use Vector temporarily to wrap the feature because that's what addSupportingFeatures takes
  Vector *tmpVec = Vector_new();
  Vector_addElement(tmpVec, feat);
  Exon_addSupportingFeatures(paddedExon, tmpVec);
  Vector_free(tmpVec);

  return paddedExon;
}

// Take the set of rough transcripts, and merge all that overlap with each other, to simplify
// the set, and hopefully make it easier to generate longer models
/*
Vector *RefineSolexaGenes_collapseTranscripts(RefineSolexaGenes *rsg, Vector *transcripts) {
  Vector *collapsedTranscripts = Vector_new();

// Make clusters
  int i;
  for (i=0; i<; i++) {
  }
  
    
// for each cluster
  for (i=0; i<; i++) {
// AIM: Make a single new transcripts which contains the collapsed exon ranges from all the others in the cluster
//   Make a range registry and add all the exon ranges from all the transcripts in to it
    RangeRegistry *reg = RangeRegistry_new();

    int j;
    for (j=0; j<; j++) {
      Transcript *trans = Gene_getTranscriptAt(gene, 0);
      int k;
      for (k=0;k<Transcript_getExonCount(trans); k++) {
        RangeRegistry_checkAndRegister(reg, (IDType)cluster, Exon_getStart(exon), Exon_getEnd(exon), 
                                       Exon_getStart(exon), Exon_getEnd(exon), 1);
      }
    }

//   Read out the ranges and make them into exons
    Vector *ranges = RangeRegistry_getRanges(reg, cluster);
    for (j=0; j<Vector_getNumElement(ranges); j++) {
      CoordPair *range = Vector_getElementAt(ranges, j);
      Exon *newEx = Exon_new();
      Exon_setStart(newEx, 
//   Add exons into new transcript

//   Add support??
//   Add into new Transcripts vector
    RangeRegistry_free(reg);
// Done cluster
  }

// Return new transcripts vector
  return collapsedTranscripts;
}
*/
 
//##################################################################
//# Containers

void RefineSolexaGenes_setRecursiveLimit(RefineSolexaGenes *rsg, int limit) {
  rsg->recursiveLimit = limit;
}

int RefineSolexaGenes_getRecursiveLimit(RefineSolexaGenes *rsg) {
  return rsg->recursiveLimit;
}

/* Not used
sub repeat_feature_adaptor {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_rfa} = $val;
  }

  return $self->{_rfa};
}
*/

SliceAdaptor *RefineSolexaGenes_getGeneSliceAdaptor(RefineSolexaGenes *rsg) {
  return rsg->geneSliceAdaptor;
}

void RefineSolexaGenes_setGeneSliceAdaptor(RefineSolexaGenes *rsg, SliceAdaptor *sa) {
  rsg->geneSliceAdaptor = sa;
}

SliceAdaptor *RefineSolexaGenes_getIntronSliceAdaptor(RefineSolexaGenes *rsg) {
  return rsg->intronSliceAdaptor;
}

void RefineSolexaGenes_setIntronSliceAdaptor(RefineSolexaGenes *rsg, SliceAdaptor *sa) {
  rsg->intronSliceAdaptor = sa;
}

/* Unused - well was used, but in a method which was unused
sub small_intron_slice_adaptor {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_sisa} = $val;
  }

  return $self->{_sisa};
}
*/

/* Not used
sub feature_cash {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_feature_cash} = $val;
  }

  return $self->{_feature_cash};
}
*/

Vector *RefineSolexaGenes_getPrelimGenes(RefineSolexaGenes *rsg) {
  return rsg->prelimGenes;
}

void RefineSolexaGenes_setPrelimGenes(RefineSolexaGenes *rsg, Vector *genes) {
  rsg->prelimGenes = genes;
}

Slice *RefineSolexaGenes_getChrSlice(RefineSolexaGenes *rsg) {
  return rsg->chrSlice;
}

void RefineSolexaGenes_setChrSlice(RefineSolexaGenes *rsg, Slice *slice) {
  rsg->chrSlice = slice;
}

Vector *RefineSolexaGenes_getIntronFeatures(RefineSolexaGenes *rsg) {
  return rsg->intronFeatures;
}

// The perl implemenetation (see below) allowed appending to the array, but I don't think that was used so I've not implemented it
void RefineSolexaGenes_setIntronFeatures(RefineSolexaGenes *rsg, Vector *features) {
  int verbosity = RefineSolexaGenes_getVerbosity(rsg);
  if (rsg->intronFeatures != NULL) {
    if (verbosity > 0) fprintf(stderr, "Trying to set intronFeatures when its already set - may want append behaviour which isn't implemented - assuming you want replace, so freeing old set\n");
    Vector_setFreeFunc(rsg->intronFeatures, DNAAlignFeature_freeImpl);
    Vector_free(rsg->intronFeatures);
  }
// Perl did this sort - it did startComp, I'm doing startEndComp
//  Vector_sort(features, SeqFeature_startCompFunc);
  Vector_sort(features, SeqFeature_startEndCompFunc);

  long maxLength = 0;
  int i;
  for (i=0;i<Vector_getNumElement(features);i++) {
    DNAAlignFeature *daf = Vector_getElementAt(features, i);
    if (DNAAlignFeature_getLength(daf) > maxLength) {
      maxLength = DNAAlignFeature_getLength(daf);
    }
  }
  if (verbosity > 0) fprintf(stderr,"Longest intron length = %ld\n",maxLength);
  RefineSolexaGenes_setLongestIntronLength(rsg, maxLength);
  
  rsg->intronFeatures = features;
}

void RefineSolexaGenes_setLongestIntronLength(RefineSolexaGenes *rsg, long maxLength) {
  rsg->longestIntronLength = maxLength;
}

long RefineSolexaGenes_getLongestIntronLength(RefineSolexaGenes *rsg) {
  return rsg->longestIntronLength;
}
/*
sub intron_features {
  my ($self, $val) = @_;

  if (defined $val) {
    push @{$self->{_introns}}, @$val;
    # make sure it is still sorted
    @{$self->{_introns}} = sort { $a->start <=> $b->start } @{$self->{_introns}};
  }
  return $self->{_introns};
}
*/

StringHash *RefineSolexaGenes_getExtraExons(RefineSolexaGenes *rsg) {
  if (rsg->extraExons == NULL) {
    rsg->extraExons = StringHash_new(STRINGHASH_MEDIUM);
  }
  return rsg->extraExons;
}

void RefineSolexaGenes_setExtraExons(RefineSolexaGenes *rsg, StringHash *extraExons) {
  rsg->extraExons = extraExons;
  rsg->extraExonsKeys = NULL;
}

char **RefineSolexaGenes_getExtraExonsKeys(RefineSolexaGenes *rsg) {
  if (rsg->extraExonsKeys == NULL) {
    rsg->extraExonsKeys  = StringHash_getKeys(RefineSolexaGenes_getExtraExons(rsg));
  }

  return rsg->extraExonsKeys;
}

ExtraExonData **RefineSolexaGenes_getExtraExonsValues(RefineSolexaGenes *rsg) {
  if (rsg->extraExonsValues == NULL) {
    rsg->extraExonsValues  = StringHash_getValues(RefineSolexaGenes_getExtraExons(rsg));
    if (rsg->extraExonsValues) 
      qsort(rsg->extraExonsValues,
            StringHash_getNumValues(RefineSolexaGenes_getExtraExons(rsg)),
            sizeof(ExtraExonData *),
            ExtraExonData_startCompFunc);
  }

  return rsg->extraExonsValues;
}




//####################################
//# config variable holders
//####################################

void RefineSolexaGenes_setIntronDb(RefineSolexaGenes *rsg, char *intronDb) {
  if (intronDb[0] != '\0') {
    if (rsg->intronDb != NULL) {
      free(rsg->intronDb);
    }
    rsg->intronDb = StrUtil_copyString(&rsg->intronDb, intronDb, 0);
  }
}

char *RefineSolexaGenes_getIntronDb(RefineSolexaGenes *rsg) {
  return rsg->intronDb;
}

void RefineSolexaGenes_setOutputDb(RefineSolexaGenes *rsg, char *outputDb) {
  if (outputDb[0] != '\0') {
    if (rsg->outputDb != NULL) {
      free(rsg->outputDb);
    }
    rsg->outputDb = StrUtil_copyString(&rsg->outputDb, outputDb, 0);
  }
}

char *RefineSolexaGenes_getOutputDb(RefineSolexaGenes *rsg) {
  return rsg->outputDb;
}

void RefineSolexaGenes_setModelDb(RefineSolexaGenes *rsg, char *modelDb) {
  if (modelDb[0] != '\0') {
    if (rsg->modelDb != NULL) {
      free(rsg->modelDb);
    }
    rsg->modelDb = StrUtil_copyString(&rsg->modelDb, modelDb, 0);
  }
}

char *RefineSolexaGenes_getModelDb(RefineSolexaGenes *rsg) {
  return rsg->modelDb;
}

void RefineSolexaGenes_setLogicNames(RefineSolexaGenes *rsg, Vector *logicNames) {
  rsg->logicNames = logicNames;
}

Vector *RefineSolexaGenes_getLogicNames(RefineSolexaGenes *rsg) {
  return rsg->logicNames;
}

void RefineSolexaGenes_setConsLims(RefineSolexaGenes *rsg, Vector *consLims) {
  rsg->consLims = consLims;
}

Vector *RefineSolexaGenes_getConsLims(RefineSolexaGenes *rsg) {
  return rsg->consLims;
}

void RefineSolexaGenes_setNonConsLims(RefineSolexaGenes *rsg, Vector *nonConsLims) {
  rsg->nonConsLims = nonConsLims;
}

Vector *RefineSolexaGenes_getNonConsLims(RefineSolexaGenes *rsg) {
  return rsg->nonConsLims;
}

void RefineSolexaGenes_setRestartNonConsLim(RefineSolexaGenes *rsg, double restartNonConsLim) {
  rsg->restartNonConsLim = restartNonConsLim;
}

double RefineSolexaGenes_getRestartNonConsLim(RefineSolexaGenes *rsg) {
  return rsg->restartNonConsLim;

}
void RefineSolexaGenes_setRetainedIntronPenalty(RefineSolexaGenes *rsg, double retainedIntronPenalty) {
  rsg->retainedIntronPenalty = retainedIntronPenalty;
}

double RefineSolexaGenes_getRetainedIntronPenalty(RefineSolexaGenes *rsg) {
  return rsg->retainedIntronPenalty;
}

void RefineSolexaGenes_setMinIntronSize(RefineSolexaGenes *rsg, int minIntronSize) {
  rsg->minIntronSize = minIntronSize;
}

int RefineSolexaGenes_getMinIntronSize(RefineSolexaGenes *rsg) {
  return rsg->minIntronSize;
}

void RefineSolexaGenes_setMaxIntronSize(RefineSolexaGenes *rsg, int maxIntronSize) {
  rsg->maxIntronSize = maxIntronSize;
}

int RefineSolexaGenes_getMaxIntronSize(RefineSolexaGenes *rsg) {
  return rsg->maxIntronSize;
}

void RefineSolexaGenes_setBestScoreType(RefineSolexaGenes *rsg, char *bestScoreType) {
  if (bestScoreType[0] != '\0') {
    if (rsg->bestScoreType != NULL) {
      free(rsg->bestScoreType);
    }
    rsg->bestScoreType = StrUtil_copyString(&rsg->bestScoreType, bestScoreType, 0);
  }
}

char *RefineSolexaGenes_getBestScoreType(RefineSolexaGenes *rsg) {
  return rsg->bestScoreType;
}

void RefineSolexaGenes_setOtherNum(RefineSolexaGenes *rsg, int otherNum) {
  rsg->otherNum = otherNum;
}

int RefineSolexaGenes_getOtherNum(RefineSolexaGenes *rsg) {
  return rsg->otherNum;
}

void RefineSolexaGenes_setOtherIsoformsType(RefineSolexaGenes *rsg, char *otherIsoformsType) {
    if (rsg->otherIsoformsType != NULL) {
      free(rsg->otherIsoformsType);
    }
  rsg->otherIsoformsType = StrUtil_copyString(&rsg->otherIsoformsType, otherIsoformsType, 0);
}

char *RefineSolexaGenes_getOtherIsoformsType(RefineSolexaGenes *rsg) {
  return rsg->otherIsoformsType;
}

void RefineSolexaGenes_setModelLogicName(RefineSolexaGenes *rsg, char *modelLN) {
  if (modelLN[0] != '\0') {
    if (rsg->modelLogicName != NULL) {
      free(rsg->modelLogicName);
    }
    rsg->modelLogicName = StrUtil_copyString(&rsg->modelLogicName, modelLN, 0);
  }
}

char *RefineSolexaGenes_getModelLogicName(RefineSolexaGenes *rsg) {
  return rsg->modelLogicName;
}

void RefineSolexaGenes_setBadModelsType(RefineSolexaGenes *rsg, char *badModelsType) {
  if (badModelsType[0] != '\0') {
    if (rsg->badModelsType != NULL) {
      free(rsg->badModelsType);
    }
    rsg->badModelsType = StrUtil_copyString(&rsg->badModelsType, badModelsType, 0);
  }
}

char *RefineSolexaGenes_getBadModelsType(RefineSolexaGenes *rsg) {
  return rsg->badModelsType;
}

void RefineSolexaGenes_setMaxNum(RefineSolexaGenes *rsg, int maxNum) {
  rsg->maxNum = maxNum;
}

int RefineSolexaGenes_getMaxNum(RefineSolexaGenes *rsg) {
  return rsg->maxNum;
}

void RefineSolexaGenes_setMaxRecursions(RefineSolexaGenes *rsg, int maxRecursions) {
  rsg->maxRecursions = maxRecursions;
}

int RefineSolexaGenes_getMaxRecursions(RefineSolexaGenes *rsg) {
  return rsg->maxRecursions;
}

void RefineSolexaGenes_setMinSingleExonLength(RefineSolexaGenes *rsg, int minSingleExonLength) {
  rsg->minSingleExonLength = minSingleExonLength;
}

int RefineSolexaGenes_getMinSingleExonLength(RefineSolexaGenes *rsg) {
  return rsg->minSingleExonLength;
}

void RefineSolexaGenes_setMinSingleExonCDSPercLength(RefineSolexaGenes *rsg, double minSingleExonCDSPercLength) {
  rsg->minSingleExonCDSPercLength = minSingleExonCDSPercLength;
}

double RefineSolexaGenes_getMinSingleExonCDSPercLength(RefineSolexaGenes *rsg) {
  return rsg->minSingleExonCDSPercLength;
}

void RefineSolexaGenes_setSingleExonModelType(RefineSolexaGenes *rsg, char *singleExonModelType) {
  if (singleExonModelType[0] != '\0') {
    if (rsg->singleExonModelType != NULL) {
      free(rsg->singleExonModelType);
    }
    rsg->singleExonModelType = StrUtil_copyString(&rsg->singleExonModelType, singleExonModelType, 0);
  }
}

char *RefineSolexaGenes_getSingleExonModelType(RefineSolexaGenes *rsg) {
  return rsg->singleExonModelType;
}

void RefineSolexaGenes_setStrictInternalSpliceSites(RefineSolexaGenes *rsg, int strictInternalSpliceSites) {
  rsg->strictInternalSpliceSites = strictInternalSpliceSites;
}

int RefineSolexaGenes_strictInternalSpliceSites(RefineSolexaGenes *rsg) {
  return rsg->strictInternalSpliceSites;
}

void RefineSolexaGenes_setStrictInternalEndSpliceSites(RefineSolexaGenes *rsg, int strictInternalEndSpliceSites) {
  rsg->strictInternalEndSpliceSites = strictInternalEndSpliceSites;
}

int RefineSolexaGenes_strictInternalEndSpliceSites(RefineSolexaGenes *rsg) {
  return rsg->strictInternalEndSpliceSites;
}

void RefineSolexaGenes_setIntronBamFiles(RefineSolexaGenes *rsg, Vector *intronBamFiles) {
  rsg->intronBamFiles = intronBamFiles;
}

Vector *RefineSolexaGenes_getIntronBamFiles(RefineSolexaGenes *rsg) {
  return rsg->intronBamFiles;
}

void RefineSolexaGenes_setAlignmentBamFile(RefineSolexaGenes *rsg, char *alignmentBamFile) {
  if (alignmentBamFile[0] != '\0') {
    if (rsg->alignmentBamFile != NULL) {
      free(rsg->alignmentBamFile);
    }
    rsg->alignmentBamFile = StrUtil_copyString(&rsg->alignmentBamFile, alignmentBamFile, 0);
  }
}

char *RefineSolexaGenes_getAlignmentBamFile(RefineSolexaGenes *rsg) {
  return rsg->alignmentBamFile;
}

void RefineSolexaGenes_setWriteIntrons(RefineSolexaGenes *rsg, int writeIntrons) {
  rsg->writeIntrons = writeIntrons;
}

int RefineSolexaGenes_writeIntrons(RefineSolexaGenes *rsg) {
  return rsg->writeIntrons;
}

void RefineSolexaGenes_setTrimUTR(RefineSolexaGenes *rsg, int trimUtr) {
  rsg->trimUtr = trimUtr;
}

int RefineSolexaGenes_trimUTR(RefineSolexaGenes *rsg) {
  return rsg->trimUtr;
}

void RefineSolexaGenes_setDryRun(RefineSolexaGenes *rsg, int dryRun) {
  rsg->dryRun = dryRun;
}

int RefineSolexaGenes_isDryRun(RefineSolexaGenes *rsg) {
  return rsg->dryRun;
}

void RefineSolexaGenes_setMax3PrimeExons(RefineSolexaGenes *rsg, int max3PrimeExons) {
  rsg->max3PrimeExons = max3PrimeExons;
}

int RefineSolexaGenes_getMax3PrimeExons(RefineSolexaGenes *rsg) {
  return rsg->max3PrimeExons;
}

void RefineSolexaGenes_setMax3PrimeLength(RefineSolexaGenes *rsg, int max3PrimeLength) {
  rsg->max3PrimeLength = max3PrimeLength;
}

int RefineSolexaGenes_getMax3PrimeLength(RefineSolexaGenes *rsg) {
  return rsg->max3PrimeLength;
}

void RefineSolexaGenes_setMax5PrimeExons(RefineSolexaGenes *rsg, int max5PrimeExons) {
  rsg->max5PrimeExons = max5PrimeExons;
}

int RefineSolexaGenes_getMax5PrimeExons(RefineSolexaGenes *rsg) {
  return rsg->max5PrimeExons;
}

void RefineSolexaGenes_setMax5PrimeLength(RefineSolexaGenes *rsg, int max5PrimeLength) {
  rsg->max5PrimeLength = max5PrimeLength;
}

int RefineSolexaGenes_getMax5PrimeLength(RefineSolexaGenes *rsg) {
  return rsg->max5PrimeLength;
}

void RefineSolexaGenes_setFilterOnOverlapThreshold(RefineSolexaGenes *rsg, int filterOnOverlapThreshold) {
  rsg->filterOnOverlapThreshold = filterOnOverlapThreshold;
}

int RefineSolexaGenes_getFilterOnOverlapThreshold(RefineSolexaGenes *rsg) {
  return rsg->filterOnOverlapThreshold;
}

void RefineSolexaGenes_setIsOneThreshold(RefineSolexaGenes *rsg, int isOneThreshold) {
  rsg->isOneThreshold = isOneThreshold;
}

int RefineSolexaGenes_getIsOneThreshold(RefineSolexaGenes *rsg) {
  return rsg->isOneThreshold;
}

void RefineSolexaGenes_setVerbosity(RefineSolexaGenes *rsg, int verbosity) {
  rsg->verbosity = verbosity;
}

int RefineSolexaGenes_getVerbosity(RefineSolexaGenes *rsg) {
  return rsg->verbosity;
}

void RefineSolexaGenes_setThreads(RefineSolexaGenes *rsg, int threads) {
  rsg->threads = threads;
}

int RefineSolexaGenes_getThreads(RefineSolexaGenes *rsg) {
  return rsg->threads;
}

void RefineSolexaGenes_setUcscNaming(RefineSolexaGenes *rsg, int ucsc_naming) {
  rsg->ucsc_naming = ucsc_naming;
}

int RefineSolexaGenes_getUcscNaming(RefineSolexaGenes *rsg) {
  return rsg->ucsc_naming;
}

void RefineSolexaGenes_setRejectIntronCutoff(RefineSolexaGenes *rsg, double rejectIntronCutoff) {
  rsg->rejectIntronCutoff = rejectIntronCutoff;
}

double RefineSolexaGenes_getRejectIntronCutoff(RefineSolexaGenes *rsg) {
  return rsg->rejectIntronCutoff;
}

void RefineSolexaGenes_setConsLim(RefineSolexaGenes *rsg, double consLim) {
  rsg->consLim = consLim;
}

double RefineSolexaGenes_getConsLim(RefineSolexaGenes *rsg) {
  return rsg->consLim;
}

void RefineSolexaGenes_setNonConsLim(RefineSolexaGenes *rsg, double nonConsLim) {
  rsg->nonConsLim = nonConsLim;
}

double RefineSolexaGenes_getNonConsLim(RefineSolexaGenes *rsg) {
  return rsg->nonConsLim;
}

void RefineSolexaGenes_setDb(RefineSolexaGenes *rsg, DBAdaptor *db) {
  rsg->db = db;
}

DBAdaptor *RefineSolexaGenes_getDb(RefineSolexaGenes *rsg) {
  return rsg->db;
}

void RefineSolexaGenes_setInputId(RefineSolexaGenes *rsg, char *inputId) {
    if (rsg->inputId != NULL) {
      free(rsg->inputId);
    }
  rsg->inputId = StrUtil_copyString(&rsg->inputId, inputId, 0);
}

char *RefineSolexaGenes_getInputId(RefineSolexaGenes *rsg) {
  return rsg->inputId;
}

void RefineSolexaGenes_setTypePrefix(RefineSolexaGenes *rsg, char *typePrefix) {
    if (rsg->typePrefix != NULL) {
      free(rsg->typePrefix);
    }
  rsg->typePrefix = StrUtil_copyString(&rsg->typePrefix, typePrefix, 0);
}

char *RefineSolexaGenes_getTypePrefix(RefineSolexaGenes *rsg) {
  return rsg->typePrefix;
}

Exon *ExonUtils_cloneExon(Exon *exon) {
  Vector *supportingFeatures = NULL;

  nExonClone++;

  Vector *origSupport = Exon_getAllSupportingFeatures(exon);
  if (origSupport != NULL && Vector_getNumElement(origSupport)) {
    supportingFeatures = Vector_new();
    int i;
    for (i=0; i<Vector_getNumElement(origSupport); i++) {
      BaseAlignFeature *baf = Vector_getElementAt(origSupport, i);

      BaseAlignFeature *newBaf = EvidenceUtils_cloneEvidence(baf);
      //fprintf(stderr,"Cloned baf %p has ref count of %d\n",newBaf, Object_getRefCount(newBaf));
      //Object_incRefCount(newBaf);
     
      Vector_addElement(supportingFeatures, newBaf);
    }
  }

  Exon *newExon = ExonUtils_createExon(Exon_getStart(exon), 
                                       Exon_getEnd(exon),
                                       Exon_getPhase(exon),
                                       Exon_getEndPhase(exon),
                                       Exon_getStrand(exon),
                                       Exon_getAnalysis(exon),
                                       supportingFeatures,
                                       Exon_getDbID(exon),
                                       Exon_getSlice(exon),
                                       Exon_getStableId(exon),
                                       Exon_getVersion(exon));

  if (supportingFeatures) Vector_free(supportingFeatures);

//  if (exon->seqCacheString) {
//    StrUtil_copyString(&newExon->seqCacheString,exon->seqCacheString,0);
//  }

  return newExon;
}


// Makes a deep copy including dbid
BaseAlignFeature *EvidenceUtils_cloneEvidence(BaseAlignFeature *feature) {
  if (feature->objectType == CLASS_DNADNAALIGNFEATURE) {
    DNAAlignFeature *newFeature = DNAAlignFeature_deepCopy(feature);
    return newFeature;
  } else if (feature->objectType == CLASS_DNAPEPALIGNFEATURE) {
    DNAPepAlignFeature *newFeature = DNAPepAlignFeature_deepCopy(feature);
    return newFeature;
  } else {
    fprintf(stderr, "ExidenceUtils:clone_Evidence Don't know what to do with a feature of type %s\n", Class_findByType(feature->objectType)->name);
    exit(EX_SOFTWARE);
  }
}

Exon *ExonUtils_createExon(long start, long end, int phase, int endPhase, int strand, Analysis *analysis, Vector *supportingFeatures, IDType dbId, Slice *slice, char *stableId, int version) {

  Exon *newExon = Exon_new();

  Exon_setStart(newExon, start);
  Exon_setEnd(newExon, end);
  Exon_setPhase(newExon, phase);
  Exon_setEndPhase(newExon, endPhase);
  Exon_setStrand(newExon, strand);
  Exon_setAnalysis(newExon, analysis);
  Exon_setDbID(newExon, dbId);
  Exon_setSlice(newExon, slice);
  if (stableId) Exon_setStableId(newExon, stableId);
  if (version)  Exon_setVersion(newExon, version);

  if (supportingFeatures) {
// Perl did this one at a time for some reason
    Exon_addSupportingFeatures(newExon, supportingFeatures);
  }

  return newExon;
}

ORFRange *ORFRange_new(long length, long start, long end) {
  ORFRange *orfRange;

  if ((orfRange = (ORFRange *)calloc(1,sizeof(ORFRange))) == NULL) {
    fprintf(stderr,"Failed allocating ORFRange\n");
    exit(EXIT_MEMORY);
  }
  orfRange->length = length;
  orfRange->start  = start;
  orfRange->end    = end;

  return orfRange;
}

int ORFRange_reverseLengthCompFunc(const void *a, const void *b) {
  ORFRange *orfRange1 = *((ORFRange **)a);
  ORFRange *orfRange2 = *((ORFRange **)b);

  return orfRange2->length - orfRange1->length;
}

void ORFRange_free(ORFRange *orfRange) {
  free(orfRange);
}

// Attempt at version of generateORFRanges which does nonMet and Met in one go
Vector *TranslationUtils_generateORFRanges(Transcript *transcript, int requireMet, int minLength, int allowReverse) {
  long len;
  int frame;
  long start;
  long end;
  char *orf;
  char *metOrf;
  int i;
  Vector *metOrfRanges = Vector_new();
  Vector *otherOrfRanges = Vector_new();

  int   lengths[6];
  char *aaSeq[6];
  char *endAaSeq[6];

  char *mRNA = Transcript_getSplicedSeq(transcript);

  int lenmRNA = strlen(mRNA);
  
  for (i=0;i<6;i++) {
    aaSeq[i] = (char *)malloc(lenmRNA/3 + 2);
  }

  long codonTableId = 1;

  //fprintf(stderr, "translateable seq = %s\n",mRNA);
  translate(mRNA, aaSeq, lengths, codonTableId, lenmRNA);

  for (i=0;i<6;i++) {
    endAaSeq[i] = aaSeq[i] + lengths[i];
    //fprintf(stderr,"aaSeq[%d] = %s\n",i,aaSeq[i]);
  }

  free(mRNA);

  int bestMetLen = 0;
  int bestOtherLen = 0;

  ORFRange *bestMetORF = NULL;
  ORFRange *bestOtherORF = NULL;

  for (frame = 0; frame < 6; frame++) {
    if (!allowReverse && frame > 2) continue;

    orf = strtok(aaSeq[frame], "*");
    while (orf != NULL && *orf != '\0') {
      metOrf = orf;
      while (*metOrf != 'M' && *metOrf != '\0') metOrf++;

      if (*orf != '\0') {
        // If have an other ORF which doesn't start at the same place as the met ORF then store that one separately to the met on
        len = -1;
        if (metOrf != orf) {
          len = strlen(orf);
          if (len > minLength && len > bestOtherLen) {
            bestOtherLen = len;
            start = (orf - aaSeq[frame]) * 3 + 1;
            if (frame < 3) {
              start += frame;
            } else {
              start -= frame-3;
            }
  
            if (frame < 3) {
              end = start + len * 3 - 1;
            } else {
              start = -1 * (start - lenmRNA - 1);
              end = start - len * 3 + 1;
            }
            fprintf(stderr, "S: %ld E: %ld T: %s\n", start, end, orf);
  
            if (orf + len != endAaSeq[frame]) { // happens when the strtok found a '*' - add 3 to end to include it in translation range
              //Vector_addElement(otherOrfRanges, ORFRange_new(len, start, end+3));
              if (!bestOtherORF) {
                bestOtherORF = ORFRange_new(len, start, end+3);
              } else {
                bestOtherORF->length = len;
                bestOtherORF->start  = start;
                bestOtherORF->end    = end+3;
              }
              //fprintf(stderr,"orf range (* ended) frame %d length %ld start %ld end %ld lenseq = %d\n", frame, len, start, end+3, lengths[frame]);
            } else {  // should be no '*' (at end of aaSeq string - this should be the last token)
              //Vector_addElement(otherOrfRanges, ORFRange_new(len, start, end));
              if (!bestOtherORF) {
                bestOtherORF = ORFRange_new(len, start, end);
              } else {
                bestOtherORF->length = len;
                bestOtherORF->start  = start;
                bestOtherORF->end    = end;
              }
              //fprintf(stderr,"orf range (NOT * ended) frame %d length %ld start %ld end %ld lenseq = %d\n", frame, len, start, end, lengths[frame]);
            }
          }
        }

        if (*metOrf != '\0') {
/* Doesn't seem to speed it up much
          if (len == -1) { 
            len = strlen(metOrf);
          } else {
            len = len - (metOrf-orf);
          }
*/
          len = strlen(metOrf);
          if (len > minLength && len > bestMetLen) {
            bestMetLen = len;
            start = (metOrf - aaSeq[frame]) * 3 + 1;
            if (frame < 3) {
              start += frame;
            } else {
              start -= frame-3;
            }
  
            if (frame < 3) {
              end = start + len * 3 - 1;
            } else {
              start = -1 * (start - lenmRNA - 1);
              end = start - len * 3 + 1;
            }
  
            if (metOrf + len != endAaSeq[frame]) { // happens when the strtok found a '*' - add 3 to end to include it in translation range
              //Vector_addElement(metOrfRanges, ORFRange_new(len, start, end+3));
              if (!bestMetORF) {
                bestMetORF = ORFRange_new(len, start, end+3);
              } else {
                bestMetORF->length = len;
                bestMetORF->start  = start;
                bestMetORF->end    = end+3;
              }
              //fprintf(stderr,"orf range (* ended) frame %d length %ld start %ld end %ld lenseq = %d\n", frame, len, start, end+3, lengths[frame]);
            } else {  // should be no '*' (at end of aaSeq string - this should be the last token)
              //Vector_addElement(metOrfRanges, ORFRange_new(len, start, end));
              if (!bestMetORF) {
                bestMetORF = ORFRange_new(len, start, end);
              } else {
                bestMetORF->length = len;
                bestMetORF->start  = start;
                bestMetORF->end    = end;
              }
              //fprintf(stderr,"orf range (NOT * ended) frame %d length %ld start %ld end %ld lenseq = %d\n", frame, len, start, end, lengths[frame]);
            }
          }
        }
      }
      orf = strtok(NULL, "*");
    }
  }


// Free ORF strings
  for (i=0; i<6; i++) {
    free(aaSeq[i]);
  }

// Sort by length
//  Vector_sort(metOrfRanges, ORFRange_reverseLengthCompFunc);
//  Vector_sort(otherOrfRanges, ORFRange_reverseLengthCompFunc);
  if (bestMetORF) fprintf(stderr, "BM S: %ld E: %ld\n", bestOtherORF->start, bestOtherORF->end);
  if (bestOtherORF) fprintf(stderr, "BO S: %ld E: %ld\n", bestOtherORF->start, bestOtherORF->end);
  if (bestMetORF) Vector_addElement(metOrfRanges, bestMetORF);
  if (bestOtherORF) Vector_addElement(otherOrfRanges, bestOtherORF);

  //fprintf(stderr,"Sorted orf ranges:\n");
  //for (i=0;i<Vector_getNumElement(orfRanges); i++) {
  //  ORFRange *orf = Vector_getElementAt(orfRanges, i);
  //  fprintf(stderr," %ld %ld %ld\n",orf->length, orf->start, orf->end);
  //}

  Vector *orfRanges = Vector_new();
  Vector_addElement(orfRanges, metOrfRanges);
  Vector_addElement(orfRanges, otherOrfRanges);

  return orfRanges;
}

/* Original
Vector *TranslationUtils_generateORFRanges(Transcript *transcript, int requireMet, int minLength, int allowReverse) {
  long len;
  int frame;
  long start;
  long end;
  char *orf;
  int i;
  Vector *orfRanges = Vector_new();

  int   lengths[6];
  char *aaSeq[6];
  char *endAaSeq[6];

  char *mRNA = Transcript_getSplicedSeq(transcript);

  int lenmRNA = strlen(mRNA);
  
  for (i=0;i<6;i++) {
    aaSeq[i] = (char *)malloc(lenmRNA/3 + 2);
  }

  long codonTableId = 1;

  //fprintf(stderr, "translateable seq = %s\n",mRNA);
  translate(mRNA, aaSeq, lengths, codonTableId, lenmRNA);

  for (i=0;i<6;i++) {
    endAaSeq[i] = aaSeq[i] + lengths[i];
    //fprintf(stderr,"aaSeq[%d] = %s\n",i,aaSeq[i]);
  }

  free(mRNA);

  for (frame = 0; frame < 6; frame++) {
    if (!allowReverse && frame > 2) continue;

    orf = strtok(aaSeq[frame], "*");
    while (orf != NULL && *orf != '\0') {
      if (requireMet) {
        while (*orf != 'M' && *orf != '\0') orf++;
      }

      if (*orf != '\0') {
        len = strlen(orf);
        if (len > minLength) {
          start = (orf - aaSeq[frame]) * 3 + 1;
          if (frame < 3) {
            start += frame;
          } else {
            start -= frame-3;
          }

          if (frame < 3) {
            end = start + len * 3 - 1;
          } else {
            //start = -1 * (start - sqinfo.len - 1);
            start = -1 * (start - lenmRNA - 1);
            end = start - len * 3 + 1;
          }

          if (orf + len != endAaSeq[frame]) { // happens when the strtok found a '*' - add 3 to end to include it in translation range
            Vector_addElement(orfRanges, ORFRange_new(len, start, end+3));
            //fprintf(stderr,"orf range (* ended) frame %d length %ld start %ld end %ld lenseq = %d\n", frame, len, start, end+3, lengths[frame]);
          } else {  // should be no '*' (at end of aaSeq string - this should be the last token)
            Vector_addElement(orfRanges, ORFRange_new(len, start, end));
            //fprintf(stderr,"orf range (NOT * ended) frame %d length %ld start %ld end %ld lenseq = %d\n", frame, len, start, end, lengths[frame]);
          }
        }
      }
      orf = strtok(NULL, "*");
    }
  }

// Free ORF strings
  for (i=0;i<6;i++) {
    free(aaSeq[i]);
  }

// Sort by length
  Vector_sort(orfRanges, ORFRange_reverseLengthCompFunc);

  //fprintf(stderr,"Sorted orf ranges:\n");
  //for (i=0;i<Vector_getNumElement(orfRanges); i++) {
  //  ORFRange *orf = Vector_getElementAt(orfRanges, i);
  //  fprintf(stderr," %ld %ld %ld\n",orf->length, orf->start, orf->end);
  //}

  return orfRanges;
}
*/

Transcript *TranslationUtils_computeTranslation(Transcript *transcript) {
  //Vector *metPredictions   = TranslationUtils_generateORFRanges(transcript, 1, 20, 0);
  //Vector *noMetPredictions = TranslationUtils_generateORFRanges(transcript, 0, 20, 0);

  Vector *predictions   = TranslationUtils_generateORFRanges(transcript, 2, 20, 0);
  Vector *metPredictions = Vector_getElementAt(predictions, 0);
  Vector *noMetPredictions = Vector_getElementAt(predictions, 1);
  Vector_free(predictions);

  //Vector *noMetPredictions = TranslationUtils_generateORFRanges(transcript, 0, 20, 0);
  // choosing the best ORF
  ORFRange *orf = NULL;

  // Here we take the best prediction with a methionine unless
  // there aren't any of the best prediction without a
  // methoinine is more than twice the length
  if (Vector_getNumElement(metPredictions) && Vector_getNumElement(noMetPredictions)) {
    ORFRange *metBest   = Vector_getElementAt(metPredictions, 0);
    ORFRange *noMetBest = Vector_getElementAt(noMetPredictions, 0);

    if (noMetBest->length > (2 * metBest->length)){
      orf = noMetBest;
    } else {
      orf = metBest;
    }
  } else if (Vector_getNumElement(metPredictions)) {
    orf = Vector_getElementAt(metPredictions, 0);
  } else if (Vector_getNumElement(noMetPredictions)) {
    orf = Vector_getElementAt(noMetPredictions, 0);
  } else {
    fprintf(stderr, "Warning: transcript has no translations\n");
    return transcript;
  }

  // add ORF to transcript
  transcript = TranslationUtils_addORFToTranscript(orf, transcript);

  Vector_setFreeFunc(metPredictions, ORFRange_free);
  Vector_setFreeFunc(noMetPredictions, ORFRange_free);
  Vector_free(metPredictions);
  Vector_free(noMetPredictions);

  return transcript;
}


char *GeneBuildUtils_getId(SeqFeature *sf, char *idStr) {
  idStr[0] = '\0';
// assume annotated seq features can have stable ids
  if (Class_isDescendent(CLASS_ANNOTATEDSEQFEATURE, sf->objectType)) {
    AnnotatedSeqFeature *asf = (AnnotatedSeqFeature *)sf;
    if (AnnotatedSeqFeature_getStableId(asf) != NULL) {
      strcpy(idStr, AnnotatedSeqFeature_getStableId(asf));
    }
  }
  if (idStr[0] == '\0') {
    sprintf(idStr,IDFMTSTR,SeqFeature_getDbID(sf));
  }

  return idStr;
}

Transcript *TranslationUtils_addORFToTranscript(ORFRange *orf, Transcript *transcript) {
  long orfStart = orf->start;
  long orfEnd   = orf->end;

  Translation *translation = Translation_new();

  //fprintf(stderr, "Best orf for %p %d long start %ld end %ld\n", transcript, orf->length, orfStart, orfEnd);

  if (orfStart > orfEnd) {
    fprintf(stderr,"orfStart > orfEnd\n");
    exit(EX_SOFTWARE);
  }

  long translationStart;
  long translationEnd;
  Exon *translationStartExon;
  Exon *translationEndExon;

  int exonCount = 0;
  long pos = 1;

  int i;
    fprintf(stderr, "Transcript %s start:%ld end:%ld CRS:%ld CRE:%ld CDS:%ld CDE:%ld NE:%ld\n", 
    Transcript_getStableId(transcript), Transcript_getStart(transcript), Transcript_getEnd(transcript), Transcript_getCodingRegionStart(transcript), Transcript_getCodingRegionEnd(transcript), Transcript_getcDNACodingStart(transcript), Transcript_getcDNACodingEnd(transcript), Transcript_getExonCount(transcript));
  for (i=0; i<Transcript_getExonCount(transcript); i++) {
    Exon *exon = Transcript_getExonAt(transcript, i);
    exonCount++;
/*
    fprintf(stderr, "exon:%d exon_length:%ld pos:%ld orf_start:%ld orf_end:%ld pos+:%ld\n", 
            exonCount, Exon_getLength(exon), pos, orfStart, orfEnd, pos + Exon_getLength(exon) - 1);
*/

    if (orfStart >= pos && orfStart <= pos + Exon_getLength(exon) - 1) {
    fprintf(stderr, "S exon:%d exon_length:%ld pos:%ld orf_start:%ld orf_end:%ld pos+:%ld\n",
            exonCount, Exon_getLength(exon), pos, orfStart, orfEnd, pos + Exon_getLength(exon) - 1);
      translationStart = orfStart - pos + 1;
      translationStartExon = exon;
    }
    if (orfEnd >= pos && orfEnd <= pos + Exon_getLength(exon) - 1) {
    fprintf(stderr, "E exon:%d exon_length:%ld pos:%ld orf_start:%ld orf_end:%ld pos+:%ld\n",
            exonCount, Exon_getLength(exon), pos, orfStart, orfEnd, pos + Exon_getLength(exon) - 1);
      translationEnd     = orfEnd - pos + 1;
      translationEndExon = exon;
    }
    pos += Exon_getLength(exon);
  }

  if (!translationStart || !translationEnd || !translationStartExon || !translationEndExon) {
    char idStr[1024];
    fprintf(stderr,"problems making the translation for %s\n", GeneBuildUtils_getId((SeqFeature *)transcript, idStr));
    fprintf(stderr,"translationStart %ld translationEnd %ld translationStartExon %p translationEndExon %p\n",
            translationStart, translationEnd, translationStartExon, translationEndExon);
    return transcript;
  } else {
    if (translationStartExon == translationEndExon && translationStart > translationEnd) {
      fprintf(stderr,"Error: translation start and end in same exon but start > end\n");
      exit(EX_SOFTWARE);
    }
    Translation_setStart(translation, translationStart);
    Translation_setEnd(translation, translationEnd);

    Translation_setStartExon(translation, translationStartExon);
    Translation_setEndExon(translation, translationEndExon);

    Transcript_setTranslation(transcript, translation);
  }

  int foundStart = 0;
  int foundEnd = 0;
  int lastEndPhase;
  int firstExon = 1;
  // print "Setting phases on transcript after adding translation\n";

  for (i=0; i<Transcript_getExonCount(transcript); i++) {
    Exon *exon = Transcript_getExonAt(transcript, i);

    Exon_setPhase(exon, -1);
    Exon_setEndPhase(exon, -1);

    // print "  Have exon " . $exon->start . " " . $exon->end . "\n";
    if (Translation_getStartExon(translation) == exon) {
      if (Translation_getStart(translation) == 1 && firstExon){
        Exon_setPhase(exon, 0);
        // print "   setting start phase on it to 0 (tstart = 1 and is start_Exon)\n";
      }
      foundStart = 1;
    } else if (foundStart && !foundEnd){
      Exon_setPhase(exon, lastEndPhase);
      // print "   setting start phase on it to last_end_phase ($last_end_phase)\n";
    }

    int endPhase;
    if (exon == Translation_getStartExon(translation)) {
      endPhase = (Exon_getEnd(exon) - (Exon_getStart(exon) + Translation_getStart(translation) - 1) +1 ) %3;
      // print "   start_Exon end phase calculated ($end_phase)\n";
    } else {
      endPhase = ((Exon_getLength(exon) + Exon_getPhase(exon)) %3);
      // print "   end phase calculated ($end_phase)\n";
    }

    if ((exon == Translation_getEndExon(translation) && Exon_getLength(exon) == Translation_getEnd(translation))){
      // print "   setting end phase to $end_phase (end exon condition)\n";
      Exon_setEndPhase(exon, endPhase);
    }

    if (exon == Translation_getEndExon(translation)) {
      foundEnd = 1;
    }

    if ((foundStart && !foundEnd)) {
      // print "   setting end phase to $end_phase (found_start and not found_end condition)\n";
      Exon_setEndPhase(exon, endPhase);
    }

    lastEndPhase = Exon_getEndPhase(exon);
    firstExon = 0;
  }

  return transcript;
}

int dumpGenes(Vector *genes, int withSupport) {
  FILE *fp = stderr;
  int i;
  int failed = 0;
  for (i=0;i<Vector_getNumElement(genes) && !failed;i++) {
    Gene *g = Vector_getElementAt(genes,i);
    fprintf(fp,"Gene %s (%s) type %s logicname %s coords: %ld %ld %d\n",Gene_getStableId(g),(Gene_getDisplayXref(g) ? DBEntry_getDisplayId(Gene_getDisplayXref(g)) : ""),Gene_getBiotype(g), Analysis_getLogicName(Gene_getAnalysis(g)), Gene_getStart(g),Gene_getEnd(g),Gene_getStrand(g));

    int j;
    for (j=0;j<Gene_getTranscriptCount(g);j++) {
      Transcript *t = Gene_getTranscriptAt(g,j);
      int k;
     
      fprintf(fp," Trans %s coords: %ld %ld %d biotype: %s logicname: %s\n",
              Transcript_getStableId(t), 
              Transcript_getStart(t),
              Transcript_getEnd(t),
              Transcript_getStrand(t),
              Transcript_getBiotype(t), 
              Transcript_getAnalysis(t) ? Analysis_getLogicName(Transcript_getAnalysis(t)) : " NO ANALYSIS!");
      if (withSupport) {
        Vector *support = Transcript_getAllSupportingFeatures(t);
        if (support) {
          for (k=0; k<Vector_getNumElement(support); k++) {
            BaseAlignFeature *baf = Vector_getElementAt(support, k);
            fprintf(fp,"   support %s coords: %ld %ld %d evalue %f hCoverage %f\n", BaseAlignFeature_getHitSeqName(baf), BaseAlignFeature_getStart(baf), BaseAlignFeature_getEnd(baf), BaseAlignFeature_getStrand(baf), BaseAlignFeature_getpValue(baf), BaseAlignFeature_gethCoverage(baf));
          }
        } else {
          fprintf(fp,"   no transcript support\n");
        }
        Vector *intronSupport = Transcript_getAllIntronSupportingEvidence(t);
        if (intronSupport) {
          for (k=0; k<Vector_getNumElement(intronSupport); k++) {
            IntronSupportingEvidence *ise = Vector_getElementAt(intronSupport, k);
            fprintf(fp,"   intron support %s type %s coords: %ld %ld %d\n", 
                    IntronSupportingEvidence_getHitName(ise), 
                    Analysis_getLogicName(IntronSupportingEvidence_getAnalysis(ise)), 
                    IntronSupportingEvidence_getStart(ise), 
                    IntronSupportingEvidence_getEnd(ise), 
                    IntronSupportingEvidence_getStrand(ise));
          }
        } else {
          fprintf(fp,"   no intron support\n");
        }
      }

      for (k=0;k<Transcript_getExonCount(t);k++) {
        Exon *e = Transcript_getExonAt(t,k);
        fprintf(fp,"  exon %s (%p) coords: %ld %ld %d\n",Exon_getStableId(e), e, Exon_getStart(e), Exon_getEnd(e), Exon_getStrand(e));
        if (withSupport) {
          Vector *support = Exon_getAllSupportingFeatures(e);
          int m;
          for (m=0; m<Vector_getNumElement(support); m++) {
            BaseAlignFeature *baf = Vector_getElementAt(support, m);
            fprintf(fp,"   support %s coords: %ld %ld %d evalue %f hCoverage %f\n", BaseAlignFeature_getHitSeqName(baf), BaseAlignFeature_getStart(baf), BaseAlignFeature_getEnd(baf), BaseAlignFeature_getStrand(baf), BaseAlignFeature_getpValue(baf), BaseAlignFeature_gethCoverage(baf));
          }
        }
      }
      Translation *tln = Transcript_getTranslation(t);
      if (tln) {
 
        fprintf(fp," translation id: %s %s %d %s %d\n",Translation_getStableId(tln), 
                Exon_getStableId(Translation_getStartExon(tln)), Translation_getStart(tln),
                Exon_getStableId(Translation_getEndExon(tln)), Translation_getEnd(tln));
        char *tSeq = Transcript_translate(t);
        fprintf(fp," translation: %s\n",tSeq);
        free(tSeq);
        Vector *tlnAttribs = Translation_getAllAttributes(tln, NULL);
        if (Vector_getNumElement(tlnAttribs)) {
          fprintf(fp, " translation attributes:\n");
          int n;
          for (n=0; n<Vector_getNumElement(tlnAttribs); n++) {
            Attribute *attrib = Vector_getElementAt(tlnAttribs, n);
            fprintf(fp, "  code %s name %s desc %s value %s\n", 
                    Attribute_getCode(attrib), 
                    Attribute_getName(attrib),
                    Attribute_getDescription(attrib),
                    Attribute_getValue(attrib));
          }
        }
      }
    }
  }
  return failed;
}

// First arg should eventually be RunnableDB type I think - need a method to get the configvar -> setting function hash
void Utilities_parseConfig(RefineSolexaGenes *rsg, config_setting_t *cfgBlock, char *label, int ignoreThrow) {
  char *DEFAULT_ENTRY_KEY = "DEFAULT";

  //my ($obj, $var_hash, $label, $ignore_throw) = @_;
  if (label == NULL) {
    fprintf(stderr, "Can't parse the config hash for object if we are give no label\n");
    exit(EX_DATAERR);
  }

  StringHash *keyCheckHash = StringHash_new(STRINGHASH_SMALL);

  int count = config_setting_length(cfgBlock);
  
  int i;
  char ucKey[2048];
  for(i = 0; i < count; ++i) {
    config_setting_t *setting = config_setting_get_elem(cfgBlock, i);

    strcpy(ucKey, config_setting_name(setting));
    StrUtil_strupr(ucKey);

    if (StringHash_contains(keyCheckHash, ucKey)) {
      fprintf(stderr, "You have two entries in your config with the same name (ignoring case) for name %s\n", ucKey);
      exit(EX_DATAERR);
    }
    StringHash_add(keyCheckHash, ucKey, setting);
  }

  // replace entries in config has with upper (perl said lower but I think its upper) case versions.
/* Don't do this modifying of the config hash - instead use the keyCheckHash to return the required section (its key is upcased already
  foreach my $k (keys %check) {
    my $old_k = $check{$k};
    my $entry = $var_hash->{$old_k};
    delete $var_hash->{$old_k};

    $var_hash->{$k} = $entry;
  }
*/

  config_setting_t *defaultSection = StringHash_getValue(keyCheckHash, DEFAULT_ENTRY_KEY);

  //config_setting_t *defaultSection = config_setting_get_member(cfgBlock, DEFAULT_ENTRY_KEY);

  if (defaultSection == NULL) {
    fprintf(stderr, "You must define a %s entry in your config", DEFAULT_ENTRY_KEY);
    exit(EX_DATAERR);
  }

  // the following will fail if there are config variables that
  // do not have a corresponding method here

  count = config_setting_length(defaultSection);

  for(i = 0; i < count; ++i) {
    config_setting_t *setting = config_setting_get_elem(defaultSection, i);

    if (!StringHash_contains(rsg->funcHash, config_setting_name(setting))) {
      fprintf(stderr, "Error: no method defined in Utilities for config variable '%s'\n", config_setting_name(setting));
      exit(EX_DATAERR);
    }
  
    SetFuncData *setFuncData = StringHash_getValue(rsg->funcHash, config_setting_name(setting));

    ConfigConverter_wrapSetCall(rsg, setFuncData, setting);
  }

  //#########################################################
  // read values of config variables for this logic name into
  // instance variable, set by method
  //#########################################################
  char ucLogic[2048];
  strcpy(ucLogic, label);
  StrUtil_strupr(ucLogic);
  
  config_setting_t *analSection = StringHash_getValue(keyCheckHash, ucLogic);

  StringHash_free(keyCheckHash, NULL);

  if (analSection != NULL) {
    // entry contains more specific values for the variables
    int count = config_setting_length(analSection);
  
    int i;
    for(i = 0; i < count; ++i) {
      config_setting_t *setting = config_setting_get_elem(analSection, i);
  
      if (!StringHash_contains(rsg->funcHash, config_setting_name(setting))) {
        fprintf(stderr, "Error: no method defined in Utilities for config variable '%s'\n", config_setting_name(setting));
        exit(EX_DATAERR);
      }
    
      SetFuncData *setFuncData = StringHash_getValue(rsg->funcHash, config_setting_name(setting));
  
      ConfigConverter_wrapSetCall(rsg, setFuncData, setting);
    }

  } else {
    if (ignoreThrow == 1 ) {
      fprintf(stderr,"Warning: Your logic_name %s doesn't appear in your config file hash - using default settings\n",  ucLogic);
    } else {
      fprintf(stderr,"Error: Your logic_name %s doesn't appear in your config file hash - using default settings\n",  ucLogic);
      exit(EX_DATAERR);
    }
  }
}

char *ConfigConverter_typeCodeToString(int code) {
  switch (code) {
    case CONFIG_TYPE_INT:
      return "INT";
      break;

    case CONFIG_TYPE_INT64:
      return "INT64";
      break;

    case CONFIG_TYPE_FLOAT:
      return "FLOAT";
      break;

    case CONFIG_TYPE_STRING:
      return "STRING";
      break;

    case CONFIG_TYPE_BOOL:
      return "BOOL";
      break;

    case CONFIG_TYPE_ARRAY:
      return "ARRAY";
      break;

    case CONFIG_TYPE_LIST:
      return "LIST";
      break;

    case CONFIG_TYPE_GROUP:
      return "GROUP";
      break;

    default:
      fprintf(stderr,"Error: Unknown type in ConfigConverter_typeCodeToString - code = %d\n", code);
      exit(EX_DATAERR);
  }
}

void ConfigConverter_wrapSetCall(RefineSolexaGenes *rsg, SetFuncData *setFuncData, config_setting_t *setting) {
  switch (config_setting_type(setting)) {
    case CONFIG_TYPE_INT:
      if (setFuncData->type != CONFIG_TYPE_INT) {
        fprintf(stderr,"Error: Inconsistency between expected config arg type and actual arg type "
                       "for %s (value = %d), type is %s expected type is %s (config type code)\n", 
                config_setting_name(setting), config_setting_get_int(setting), 
                ConfigConverter_typeCodeToString(config_setting_type(setting)), 
                ConfigConverter_typeCodeToString(setFuncData->type));
        exit(EX_DATAERR);
      }
      setFuncData->setFunc.setIntValue(rsg, config_setting_get_int(setting));
      break;

    case CONFIG_TYPE_INT64:
      if (setFuncData->type != CONFIG_TYPE_INT64) {
        fprintf(stderr,"Error: Inconsistency between expected config arg type and actual arg type "
                       "for %s (value = "IDFMTSTR"), type is %s expected type is %s (config type code)\n", 
                config_setting_name(setting), config_setting_get_int64(setting),
                ConfigConverter_typeCodeToString(config_setting_type(setting)), 
                ConfigConverter_typeCodeToString(setFuncData->type));
        exit(EX_DATAERR);
      }
      setFuncData->setFunc.setInt64Value(rsg, config_setting_get_int64(setting));
      break;

    case CONFIG_TYPE_FLOAT:
      if (setFuncData->type != CONFIG_TYPE_FLOAT) {
        fprintf(stderr,"Error: Inconsistency between expected config arg type and actual arg type "
                       "for %s (value = %f), type is %s expected type is %s (config type code)\n", 
                config_setting_name(setting), config_setting_get_float(setting),
                ConfigConverter_typeCodeToString(config_setting_type(setting)), 
                ConfigConverter_typeCodeToString(setFuncData->type));
        exit(EX_DATAERR);
      }
      setFuncData->setFunc.setFloatValue(rsg, config_setting_get_float(setting));
      break;

    case CONFIG_TYPE_STRING:
      if (setFuncData->type != CONFIG_TYPE_STRING) {
        fprintf(stderr,"Error: Inconsistency between expected config arg type and actual arg type "
                       "for %s (value = %s), type is %s expected type is %s (config type code)\n", 
                config_setting_name(setting), config_setting_get_string(setting),
                ConfigConverter_typeCodeToString(config_setting_type(setting)), 
                ConfigConverter_typeCodeToString(setFuncData->type));
        exit(EX_DATAERR);
      }
      //if (config_setting_get_string(setting)[0] != '\0') {
        setFuncData->setFunc.setStringValue(rsg, config_setting_get_string(setting));
     // }
      break;

    case CONFIG_TYPE_BOOL:
      if (setFuncData->type != CONFIG_TYPE_BOOL) {
        fprintf(stderr,"Error: Inconsistency between expected config arg type and actual arg type "
                       "for %s (value = %d), type is %s expected type is %s (config type code)\n", 
                config_setting_name(setting), config_setting_get_bool(setting),
                ConfigConverter_typeCodeToString(config_setting_type(setting)), 
                ConfigConverter_typeCodeToString(setFuncData->type));
        exit(EX_DATAERR);
      }
      setFuncData->setFunc.setBoolValue(rsg, config_setting_get_bool(setting));
      break;

    case CONFIG_TYPE_ARRAY:
      if (setFuncData->type != CONFIG_TYPE_ARRAY) {
        fprintf(stderr,"Error: Inconsistency between expected config arg type and actual arg type "
                       "for %s, type is %s expected type is %s (config type code)\n", 
                config_setting_name(setting), 
                ConfigConverter_typeCodeToString(config_setting_type(setting)), 
                ConfigConverter_typeCodeToString(setFuncData->type));
        exit(EX_DATAERR);
      }
      ConfigConverter_wrapArraySetCall(rsg, setFuncData, setting);
      break;

    case CONFIG_TYPE_LIST:
      if (setFuncData->type != CONFIG_TYPE_LIST) {
        fprintf(stderr,"Error: Inconsistency between expected config arg type and actual arg type "
                       "for %s, type is %s expected type is %s (config type code)\n", 
                config_setting_name(setting), 
                ConfigConverter_typeCodeToString(config_setting_type(setting)), 
                ConfigConverter_typeCodeToString(setFuncData->type));
        exit(EX_DATAERR);
      }
      ConfigConverter_wrapListSetCall(rsg, setFuncData, setting);
      break;

    case CONFIG_TYPE_GROUP:
      if (setFuncData->type != CONFIG_TYPE_GROUP) {
        fprintf(stderr,"Error: Inconsistency between expected config arg type and actual arg type "
                       "for %s, type is %s expected type is %s (config type code)\n", 
                config_setting_name(setting), 
                ConfigConverter_typeCodeToString(config_setting_type(setting)), 
                ConfigConverter_typeCodeToString(setFuncData->type));
        exit(EX_DATAERR);
      }
      ConfigConverter_wrapGroupSetCall(rsg, setFuncData, setting);
      break;
    
    default:
      fprintf(stderr,"Error: Unknown type in config type switch for %s, type is %d expected type is %d (config type code)\n", 
              config_setting_name(setting), config_setting_type(setting), setFuncData->type);
      exit(EX_DATAERR);
  }
}

void ConfigConverter_wrapGroupSetCall(RefineSolexaGenes *rsg, SetFuncData *setFuncData, config_setting_t *setting) {
  fprintf(stderr,"NIY: Don't know how to handle Group config with name %s\n", config_setting_name(setting));
  exit(EX_DATAERR);
}

void ConfigConverter_wrapListSetCall(RefineSolexaGenes *rsg, SetFuncData *setFuncData, config_setting_t *setting) {
  if (setFuncData->subFunc) {
    setFuncData->subFunc(rsg, setting);
  } else {
    fprintf(stderr,"Don't know how to handle List config with name %s\n", config_setting_name(setting));
    exit(EX_DATAERR);
  }
}

void ConfigConverter_wrapArraySetCall(RefineSolexaGenes *rsg, SetFuncData *setFuncData, config_setting_t *setting) {
  Vector *results = Vector_new();
  if (setFuncData->subType == CONFIG_TYPE_STRING) {
    char *tmp;
    int i;
    for (i=0; i<config_setting_length(setting); i++) {
      const char *val = config_setting_get_string_elem(setting, i);
      Vector_addElement(results, StrUtil_copyString(&tmp, (char*)val, 0));
    }
  } else if (setFuncData->subType == CONFIG_TYPE_FLOAT) {
    int i;
    for (i=0; i<config_setting_length(setting); i++) {
      double val = config_setting_get_float_elem(setting, i);
      Vector_addElement(results, double_new(val));
    }
  } else {
    fprintf(stderr,"Array type other than string not implemented yet\n");
    exit(EX_DATAERR);
  }
  setFuncData->setFunc.setVectorValue(rsg, results);
}

void RefineSolexaGenes_parseIntronBamFilesConfig(RefineSolexaGenes *rsg, config_setting_t *setting) {
  Vector *intronBamFiles = Vector_new(); 

  int i;
  for (i=0; i<config_setting_length(setting); i++) {
    config_setting_t *section = config_setting_get_elem(setting, i);

    if (config_setting_type(section) != CONFIG_TYPE_GROUP) {
      fprintf(stderr,"Error: Expected a CONFIG_TYPE_GROUP setting within the IntronBamFiles block\n");
      exit(EX_DATAERR);
    }

    const char *file;
    int mixedBam;
    int depth;
    char *tmp;
    if (!(config_setting_lookup_string(section, "FILE", &file) &&
          config_setting_lookup_int(section, "MIXED_BAM", &mixedBam) &&
          config_setting_lookup_int(section, "DEPTH", &depth))) {
      fprintf(stderr,"Error: Missing at least one required arg (FILE, MIXED_BAM or DEPTH) in INTRON_BAM_FILES config block\n");
      exit(EX_DATAERR);
    }

    Vector *groupNames = NULL;
    config_setting_t *groupNameSetting = config_setting_get_member(section, "GROUPNAME");
    if (!groupNameSetting) {
      fprintf(stderr,"Error: Missing GROUPNAME required arg in INTRON_BAM_FILES config block\n");
      exit(EX_DATAERR);
    }
    
    int nGroup = config_setting_length(groupNameSetting);
    if (nGroup) {
      groupNames = Vector_new();
      
      int j;
      for (j=0; j<config_setting_length(groupNameSetting); j++) {
        const char *val = config_setting_get_string_elem(groupNameSetting, j);
        Vector_addElement(groupNames, StrUtil_copyString(&tmp, (char *)val, 0));
      }
    }
    Vector_addElement(intronBamFiles, IntronBamConfig_new((char *)file, mixedBam, depth, groupNames));
  }

  RefineSolexaGenes_setIntronBamFiles(rsg, intronBamFiles);
}


/*
=head2 select_random_db 

  Arg       : String - HashKey pointing to an entry in Databases.pm - 
              Either this key points to a key in the hash %DATABASES or in %DISTRIBUTED_DBS 

  Fuction   : The fuction reads the 2 hashes %DATABASES and %DISTRIBUTED_DBS which are exported by 
              Databases.pm. If the Argument is found in  %DATABASES, the name is returned. 
              if the Argument points to an entry in %DISTRIBUTED_DBS and %DISTRIBUTED_DBS{$arg} is 
              an array reference, an element is randomly picked out of this array and returned. 

              This function is basically used to spread the load over different db servers randomly. 
              with the $DISTRIBUTED_DBS array in Databases.pm 
 
  Returntype: String 

=cut
*/
/* NIY
char * select_random_db {  
  my ( $name ) = @_; 
 
  my $tmp;    
  if (exists $DATABASES->{$name} && ref($DATABASES->{$name}) =~m/AREF/  ) { 
    $tmp = $DATABASES ; 
  } elsif ( exists  $DISTRIBUTED_DBS->{$name} && ref($DISTRIBUTED_DBS->{$name}) =~m/ARRAY/   ) {  
    $tmp = $DISTRIBUTED_DBS; 
  } 
  if ( defined $tmp ) { 
     my @array = @{ $tmp->{$name} };   
     my $randomIndex = rand(@array);  
     $name = $array[$randomIndex]; 
     print "Random database selected : $name \n"; 
  }
  return $name ; 
}
*/

/*
=head2 get_dbadaptor

  Arg [0]   : Bio::EnsEMBL::Analysis::RunnableDB
  Arg [1]   : String - key of database hash
  Arg [2]   : return a non-standard adaptor [ valie values : 'pipeline' 'compara' 'functgenomics' or undef ] 
  Arg [3]   : flag to attch dna_db nor not 

  Function  : Returns a Bio::EnsEMBL::DBSQL::DBAdaptor for a given hash key.
              or a Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor if requested
              Requires proper configuration of
              Bio::EnsEMBL::Analysis::Config::Databases

  Returntype: Bio::EnsEMBL:DBSQL::DBAdaptor or Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
              Returns undef if an entry for the database is found but is empty.
  Exceptions: throw if key can't be found in Databases.pm

=cut
*/

// Hacky method for getting DbAdaptor given a string
// Now replaced with more complete version
#if 0
DBAdaptor *RefineSolexaGenes_getDbAdaptor(RefineSolexaGenes *rsg, char *alias) {
  if (rsg->adaptorAliasHash == NULL) {
    fprintf(stderr, "Error: adaptorAliasHash is NULL - no adaptor aliases set up - bye\n");
    exit(1);
  }
  if (StringHash_contains(rsg->adaptorAliasHash, alias)) {
// Not sure if this is what I need 
    return StringHash_getValue(rsg->adaptorAliasHash, alias);
  } else {
    fprintf(stderr,"Error: No database with alias %s\n", alias);
    exit(1);
  }
}
#endif

void RefineSolexaGenes_setDatabaseConfig(RefineSolexaGenes *rsg, config_setting_t *setting) {
  rsg->databaseConfig = setting;
}

config_setting_t *RefineSolexaGenes_getDatabaseConfig(RefineSolexaGenes *rsg) {
  return rsg->databaseConfig;
}

DBAdaptor *BaseGeneBuild_getDbAdaptor(RefineSolexaGenes *rsg, char *alias, int isNonStandard, int dontUseDnaDb) {
  DBAdaptor *db = NULL;
  int tryToAttachDnaDb = 0;
  int verbosity = RefineSolexaGenes_getVerbosity(rsg);

  config_setting_t *setting = RefineSolexaGenes_getDatabaseConfig(rsg);

  config_setting_t *databases = config_setting_get_member(setting, "DATABASES");

  const char *DNA_DBNAME;
  config_setting_lookup_string(setting, "DNA_DBNAME", &DNA_DBNAME);

  // NIY $alias = select_random_db($alias);  
  if (rsg->adaptorAliasHash == NULL) {
    fprintf(stderr, "Error: adaptorAliasHash is NULL - should have been allocated by now - bye\n");
    exit(EX_SOFTWARE);
  }

  StringHash *hash = rsg->adaptorAliasHash;

  if (!StringHash_contains(hash, alias)) { // if we don't already have an entry for this ...
    config_setting_t *section = config_setting_get_member(databases, alias);
    if (section != NULL) { 
/*
      my $constructor_args = $DATABASES->{$name};

      if ( scalar( keys( %{$constructor_args} ) ) == 0 ) {
        // The entry is empty.  Warn about this, but don't throw.
        // Return undef.
        warning(
             sprintf( "Empty entry for database '%s' in Databases.pm\n",
                      $name ) );
        return undef;
      }
*/

      // check if we got all arguments
      const char *user;
      const char *host;
      const char *dbName;
      int   port;
      if (!(config_setting_lookup_string(section, "user", &user) &&
            config_setting_lookup_string(section, "dbname", &dbName) &&
            config_setting_lookup_string(section, "host", &host) &&
            config_setting_lookup_int(section, "port", &port))) {
        fprintf(stderr,"Error: Missing at least one required arg (user, dbname, host or port) in DATABASES config block for %s\n", alias);
        exit(EX_DATAERR);
      }

/* NIY
      if ( defined $non_standard_db_adaptor ) { // value of 
        if (    $non_standard_db_adaptor =~ m/1/ || $non_standard_db_adaptor eq "pipeline" ) {
          require Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
          $db = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new( %$constructor_args );

        } elsif ( $non_standard_db_adaptor =~ m/compara/ ) {
          require Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
          unless ( $$constructor_args{'-species'} ) {
            throw("need species !\n");
          }
          $db =
            Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new( %$constructor_args );

        } elsif (    $non_standard_db_adaptor =~ m/functgenomics/i
                  || $non_standard_db_adaptor =~ m/funcgen/i )
        {
          require Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
          # funcgen adaptor needs species
          if ( $$constructor_args{'-species'} ) {
            $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
                                                        %$constructor_args, );
          } else {
            throw( "if you require a Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor "
                   . "you need to provide a -speices flag.\n" );
          }
        }
      } else {
*/
        const char *pass;
        config_setting_lookup_string(section, "pass", &pass);
        db = DBAdaptor_new((char *)host, (char *)user, (char *)pass, (char *)dbName, port, NULL);
       // $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(%$constructor_args);
        // it's a core db so try to attach dna_db
        tryToAttachDnaDb = 1;
//      }

      // this bit is attaching a dna_db .
      if (DNA_DBNAME != NULL && strcmp(alias, DNA_DBNAME) && tryToAttachDnaDb) {

        if (dontUseDnaDb) {
          if (verbosity > 0) fprintf(stderr,"\nNot attaching a DNA_DB to %s\n.", alias);
          // if two different species are considered, the
          // not_use_dna_database is set to 1 to avoid adding the second species to
          // the first one

        } else {
          // there's a little danger if we have multiple diffeent
          // species in our "Databases.pm" file. We need to avoid that the wrong
          // dna db is attached, ie a mouse core with a human dna db.

          if (verbosity > 0) fprintf(stderr,"\nAttaching DNA_DB %s to %s...\n", DNA_DBNAME, alias);
          if (DNA_DBNAME[0] == '\0') {  
            fprintf(stderr, "You're using an empty string as dna_dbname in your Databases config file\n"); 
            exit(EX_DATAERR);
          } 
          DBAdaptor *dnaDb = BaseGeneBuild_getDbAdaptor(rsg, (char *)DNA_DBNAME, 0, 0);

          // try to get default asm+ species name for OTHER db - does not work
          // for comapra database
          MetaContainer *coreMC = DBAdaptor_getMetaContainer(db);

          char *coreDbAsm     = MetaContainer_getDefaultAssembly(coreMC);
          //Species *coreDbSpecies = MetaContainer_getSpecies(coreMC);

          // get the same for dna-db
          MetaContainer *dnaMC = DBAdaptor_getMetaContainer(dnaDb);

          char *dnaDbAsm     = MetaContainer_getDefaultAssembly(dnaMC);
          //Species *dnaDbSpecies = MetaContainer_getSpecies(dnaMC);

          int dbsAreCompatible = 1;

          if (strcmp(coreDbAsm, dnaDbAsm)) { // assemblies differ
            fprintf(stderr, "You're trying to add a DNA_DB with assembly %s to "
                            "a core/cdna/otherfeatures DB with assembly %s ...\n\t"
                            "that's incompatbile. I will not add dna_database "
                            "%s to core %s\n", dnaDbAsm, coreDbAsm, DBConnection_getDbName(dnaDb->dbc), DBConnection_getDbName(db->dbc));

            dbsAreCompatible = 0;
          }

          /* If the assemblies are the same why should the species be different?
          if (strcmp(Species_getCommonName(coreDbSpecies), Species_getCommonName(dnaDbSpecies))) {  // species are different
            fprintf(stderr, "You're trying to add a DNA_DB with species %s to "
                            "a core database with speices: %s - this does not work\n"
                            "try to not use any DNA_DATABASE name in your Databases config file\n", 
                    Species_getCommonName(dnaDbSpecies), Species_getCommonName(coreDbSpecies));
            dbsAreCompatible = 0;
          }
          */
          //Species_free(dnaDbSpecies);
          //Species_free(coreDbSpecies);
          free(coreDbAsm);
          free(dnaDbAsm);

          if (dbsAreCompatible) {
            DBAdaptor_setDNADBAdaptor(db, dnaDb);
            fprintf(stderr,"\nAttaching DNA_DB %s to %s\n", DBConnection_getDbName(dnaDb->dbc), DBConnection_getDbName(db->dbc));
          }
        }
      } else {
        if ( !strcmp(alias, DNA_DBNAME)) {
          fprintf(stderr, "\nNot attaching DNA_DB to %s which has DNA_DBNAME...\n", alias); 
        } else {
          fprintf(stderr, "You haven't defined a DNA_DBNAME in your Databases config file\n");
        }
      }
    } else {
      fprintf(stderr, "No entry in Databases config file hash for %s\n", alias);
      exit(EX_DATAERR);
    }

    StringHash_add(hash, alias, db);
  } else {
    db = StringHash_getValue(hash, alias);
  }

  return db;
}

void RunnableDB_readDatabaseConfig(RefineSolexaGenes *rsg, char *configFile) {
  config_t cfg;

  config_init(&cfg);

  /* Read the file. If there is an error, report it and exit. */
  if (!config_read_file(&cfg, configFile)) {
    fprintf(stderr, "%s:%d - %s\n", config_error_file(&cfg),
            config_error_line(&cfg), config_error_text(&cfg));
    config_destroy(&cfg);
    exit(EX_NOINPUT);
  }

  config_setting_t *cfgBlock = config_lookup(&cfg, "Config");
  if (cfgBlock == NULL) {
    fprintf(stderr,"Missing config block 'Config'\n");
  }

  RefineSolexaGenes_setDatabaseConfig(rsg, cfgBlock);
}

int SeqFeat_lengthCompFunc(const void *a, const void *b) {
  SeqFeature *e1 = *((SeqFeature **)a);
  SeqFeature *e2 = *((SeqFeature **)b);

  if (SeqFeature_getLength(e1) > SeqFeature_getLength(e2)) {
    return 1;
  } else if (SeqFeature_getLength(e1) < SeqFeature_getLength(e2)) {
    return -1;
  } else {
    return 0;
  }
}
