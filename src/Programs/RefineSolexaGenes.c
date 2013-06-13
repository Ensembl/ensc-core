#include "RefineSolexaGenes.h"
#include <stdio.h>

#include "EnsC.h"

#include "DBAdaptor.h"
#include "CoordSystemAdaptor.h"
#include "AnalysisAdaptor.h"
#include "GeneAdaptor.h"
#include "DNAAlignFeatureAdaptor.h"
#include "BaseAdaptor.h"
#include "Vector.h"
#include "Slice.h"
#include "StrUtil.h"
#include "DNAAlignFeature.h"
#include "Gene.h"
#include "Intron.h"
#include "IntronSupportingEvidence.h"
#include "Transcript.h"

#include "sam.h"
#include "bam.h"
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
static int  limit = 0;

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

CigarBlock *CigarBlock_new(CigarBlockType type, long start, long end) {
  CigarBlock *block;

  if ((block = (CigarBlock *)calloc(1,sizeof(CigarBlock))) == NULL) {
    fprintf(stderr,"Failed allocating CigarBlock\n");
    exit(1);
  }
  block->type  = type;
  block->start = start;
  block->end   = end;

  return block;
}

int CigarBlock_startCompFunc(const void *a, const void *b) {
  CigarBlock *cb1 = *((CigarBlock **)a);
  CigarBlock *cb2 = *((CigarBlock **)b);

  return cb1->start - cb2->start;
}

void CigarBlock_free(CigarBlock *block) {
  free(block);
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
  SeqFeature *sf2 = Vector_getElementAt(mod1->features, 0);
  return SeqFeature_getStart(sf1) - SeqFeature_getStart(sf2);
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
    exit(1);
  }

  return ted;
}

void TranscriptExtraData_free(TranscriptExtraData *ted) {
  free(ted);
}

typedef struct ExtraExonDataStruct {
  int nCoord;
  long *coords;
  int score;
} ExtraExonData;

ExtraExonData *ExtraExonData_new(long *coords, int nCoord) {
  ExtraExonData *eed;

  if ((eed = (ExtraExonData *)calloc(1,sizeof(ExtraExonData))) == NULL) {
    fprintf(stderr,"Failed allocating ExtraExonData\n");
    exit(1);
  }

  eed->nCoord = nCoord;

  if ((eed->coords = (long *)calloc(nCoord,sizeof(long))) == NULL) {
    fprintf(stderr,"Failed allocating ExtraExonData coords for %d coords\n", nCoord);
    exit(1);
  }

  memcpy(eed->coords, coords, nCoord*sizeof(long));

  return eed;
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
    exit(1);
  }
  ic->prevExonEnd   = prevExonEnd;
  ic->nextExonStart = nextExonStart;
  ic->strand        = strand;
  ic->isCanonical   = isCanon;
  ic->score         = score;

  return ic;
}

void IntronCoords_free(IntronCoords *ic) {
  free(ic);
}

ModelCluster *ModelCluster_new() {
  ModelCluster *mc;

  if ((mc = (ModelCluster *)calloc(1,sizeof(ModelCluster))) == NULL) {
    fprintf(stderr,"Failed allocating ModelCluster\n");
    exit(1);
  }

  return mc;
}

void ModelCluster_free(ModelCluster *mc) {
  free(mc);
}

#ifdef RSG_DRIVER
int main(int argc, char *argv[]) {
  RefineSolexaGenes *rsg = RefineSolexaGenes_new(NULL);

/*
  RefineSolexaGenes_setIntronDb(rsg, char *intronDb);
  RefineSolexaGenes_setOutputDb(rsg, char *outputDb);
  RefineSolexaGenes_setModelDb(rsg, char *modelDb);
  RefineSolexaGenes_setLogicNames(rsg, Vector *logicNames);
  RefineSolexaGenes_setRetainedIntronPenalty(rsg, double retainedIntronPenalty);
  RefineSolexaGenes_setMinIntronSize(rsg, int minIntronSize);
  RefineSolexaGenes_setMaxIntronSize(rsg, int maxIntronSize);
  RefineSolexaGenes_setBestScoreType(rsg, char *bestScoreType);
  RefineSolexaGenes_setOtherNum(rsg, int otherNum);
  RefineSolexaGenes_setOtherNum(rsg, int otherNum);
  RefineSolexaGenes_setOtherIsoformsType(rsg, char *otherIsoformsType);
  RefineSolexaGenes_setModelLogicName(rsg, char *modelLN);
  RefineSolexaGenes_setBadModelsType(rsg, char *badModelsType);
  RefineSolexaGenes_setMaxNum(rsg, int maxNum);
  RefineSolexaGenes_setMaxRecursions(rsg, int maxRecursions);
  RefineSolexaGenes_setMinSingleExonLength(rsg, int minSingleExonLength);
  RefineSolexaGenes_setMinSingleExonCDSLength(rsg, int minSingleExonCDSLength);
  RefineSolexaGenes_setSingleExonModelType(rsg, char *singleExonModelType);
  RefineSolexaGenes_setStrictInternalSpliceSites(rsg, int strictInternalSpliceSites);
  RefineSolexaGenes_setStrictInternalEndSpliceSites(rsg, int strictInternalEndSpliceSites);
  RefineSolexaGenes_setIntronBamFiles(rsg, Vector *intronBamFiles);
  RefineSolexaGenes_setWriteIntrons(rsg, int writeIntrons);
  RefineSolexaGenes_setTrimUtr(rsg, int trimUtr);
  RefineSolexaGenes_setMax3PrimeExons(rsg, int max3PrimeExons);
  RefineSolexaGenes_setMax3PrimeLength(rsg, int max3PrimeLength);
  RefineSolexaGenes_setMax5PrimeExons(rsg, int max5PrimeExons);
  RefineSolexaGenes_setMax5PrimeLength(rsg, int max5PrimeLength);
  RefineSolexaGenes_setFilterOnOverlapThreshold(rsg, int filterOnOverlapThreshold);
  RefineSolexaGenes_setRejectIntronCutoff(rsg, double rejectIntronCutoff);
*/

  RefineSolexaGenes_fetchInput(rsg);
  RefineSolexaGenes_run(rsg);

//  RefineSolexaGenes_writeOutput();
}
#endif

RefineSolexaGenes *RefineSolexaGenes_new(char *configFile) {
  RefineSolexaGenes *rsg;
  
  if ((rsg = calloc(1,sizeof(RefineSolexaGenes))) == NULL) {
    fprintf(stderr,"Failed allocating RefineSolexaGenes\n");
    exit(1);
  }

// Hack for now to allow me to run without configuration 
  if (configFile) {
    fprintf(stderr, "Error: Config file specified but config reading not implemented - bye!\n");
    exit(1);
    //RunnableDB_readAndCheckConfig(rsg, configFile);
  } else {
    fprintf(stderr, "WARNING: Running without reading config (config reading not implemented)\n");
  }

//  $self->read_and_check_config($REFINESOLEXAGENES_CONFIG_BY_LOGIC);

  // Hard limit to the number of possible paths to explore
  RefineSolexaGenes_setRecursiveLimit(rsg, 10000);

  // initialise intron feature cash
  // Doesn't seem to be used
  //my %feature_hash;
  //$self->feature_cash(\%feature_hash);

  return rsg;
}

// Hacky method for getting DbAdaptor given a string
DBAdaptor *RefineSolexaGenes_getDbAdaptor(RefineSolexaGenes *rsg, char *alias) {
  if (StringHash_contains(rsg->adaptorAliasHash, alias)) {
// Not sure if this is what I need 
    return StringHash_getValue(rsg->adaptorAliasHash, alias);
  } else {
    fprintf(stderr,"Error: No database with alias %s\n", alias);
    exit(1);
  }
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
  char *logic = Analysis_getLogicName(RefineSolexaGenes_getAnalysis(rsg));
  // fetch adaptors and store them for use later
  // explicitly attach the ref db

  RefineSolexaGenes_setDb(rsg, RefineSolexaGenes_getDbAdaptor(rsg, "REFERENCE_DB"));

  DBAdaptor *db = RefineSolexaGenes_getDb(rsg);

  if (RefineSolexaGenes_getIntronDb(rsg)) {
    RefineSolexaGenes_setIntronSliceAdaptor(rsg, DBAdaptor_getSliceAdaptor(RefineSolexaGenes_getDbAdaptor(rsg, RefineSolexaGenes_getIntronDb(rsg))));
  }
// Unused  $self->repeat_feature_adaptor($self->db->get_RepeatFeatureAdaptor);

  RefineSolexaGenes_setGeneSliceAdaptor(rsg, DBAdaptor_getSliceAdaptor(RefineSolexaGenes_getDbAdaptor(rsg, RefineSolexaGenes_getModelDb(rsg))));

  // want a slice and a full chromsome to keep everything on in the same coords
  Slice *slice = RefineSolexaGenes_fetchSequence(rsg, RefineSolexaGenes_getInputId());
  Slice *geneSlice =  SliceAdaptor_fetchByRegion(RefineSolexaGenes_getGeneSliceAdaptor(rsg),
                                                 "toplevel",
                                                 Slice_getSeqRegionName(slice),
                                                 Slice_getStart(slice),
                                                 Slice_getEnd(slice),
                                                 1,
                                                 NULL,
                                                 0);

  Slice *chrSlice = SliceAdaptor_fetchByRegion(DBAdaptor_getSliceAdaptor(db), 
                                               "toplevel", 
                                               Slice_getSeqRegionName(slice),
                                               POS_UNDEF,
                                               POS_UNDEF,
                                               1,
                                               NULL,
                                               0);

  RefineSolexaGenes_setChrSlice(rsg, chrSlice);

  // we want to fetch and store all the rough transcripts at the beginning - speeds things up
  // also we want to take out tiny introns - merge them into longer structures
  Vector *genes;
  char *modelLogicName = RefineSolexaGenes_getModelLogicName(rsg);
  if (modelLogicName != NULL) {
    genes = Slice_getAllGenes(geneSlice, NULL, modelLogicName, 1, NULL, NULL);
    fprintf(stderr,"Got %d genes with logic name %s\n", Vector_getNumElement(genes), modelLogicName);
  } else {
    genes = Slice_getAllGenes(geneSlice, NULL, NULL, 1, NULL, NULL);
    fprintf(stderr, "Got %d genes\n",Vector_getNumElement(genes)); 
  }

  Vector *prelimGenes = Vector_new();
  int i;
  for (i=0; Vector_getNumElement(genes); i++) {
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
    long gc = (overlap / Gene_getLength(gene) * 1000) / 10;
    long sc =  (overlap / Slice_getLength(slice) * 1000) /10;
    printf("Gene has %ld %% overlap with the slice\nSlice has %ld %% overlap with the gene\n", gc, sc);
    if ( gc <= 10 && sc <= 10) {
      printf("Rejecting\n");
      continue;
    } 
    Vector_addElement(prelimGenes, gene);
  }
  fprintf(stderr, "Got %d genes after filtering boundary overlaps\n", Vector_getNumElement(prelimGenes)); 

  // determine strandedness ( including splitting merged genes )
  RefineSolexaGenes_setPrelimGenes(rsg, prelimGenes);

  Vector *intronBamFiles = RefineSolexaGenes_getIntronBamFiles(rsg);

  if (Vector_getNumElement(intronBamFiles) > 0) {
    int i;
    for (i=0; i<Vector_getNumElement(intronBamFiles); i++) {
      // Perl returned a hash but only seemed to use the FILE element
      IntronBamConfig *intronBamConf = Vector_getElementAt(intronBamFiles, i);
      char *intronFile = intronBamConf->fileName;
      char region[2048];
      int ref;
      int begRange;
      int endRange;
    
#ifdef _PBGZF_USE
      bam_set_num_threads_per(5);
#endif
      samfile_t *sam = samopen(intronFile, "rb", 0);
      if (sam == NULL) {
        fprintf(stderr, "Bam file %s not found\n", intronFile);
        exit(1);
      }

      bam_index_t *idx;
      idx = bam_index_load(intronFile); // load BAM index
      if (idx == 0) {
        fprintf(stderr, "BAM index file is not available.\n");
        exit(1);
      }

      long long count = 0;
      sprintf(region,"%s:%ld-%ld", Slice_getSeqRegionName(slice),
                                   Slice_getSeqRegionStart(slice),
                                   Slice_getSeqRegionEnd(slice));
      bam_parse_region(sam->header, region, &ref, &begRange, &endRange);
      if (ref < 0) {
        fprintf(stderr, "Invalid region %s %ld %ld\n", Slice_getSeqRegionName(slice),
                                                       Slice_getSeqRegionStart(slice),
                                                       Slice_getSeqRegionEnd(slice));
        exit(1);
      }

      // need to seamlessly merge here with the dna2simplefeatures code
      RefineSolexaGenes_bamToIntronFeatures(rsg, intronBamConf, sam, idx, ref, begRange, endRange);

      bam_index_destroy(idx);
      samclose(sam);
    }
  } else {
    // pre fetch all the intron features
    RefineSolexaGenes_dnaToIntronFeatures(rsg, Slice_getStart(slice), Slice_getEnd(slice));
// NIY ??
//    $self->intron_slice_adaptor->dbc->disconnect_when_inactive(1);
  }
// NIY ??
//  $self->db->disconnect_when_inactive(1);
//  $self->gene_slice_adaptor->dbc->disconnect_when_inactive(1);
}

void  RefineSolexaGenes_run(RefineSolexaGenes *rsg) {
  RefineSolexaGenes_refineGenes(rsg);
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

  int i;
// GENE:
  for (i=0; i<Vector_getNumElement(prelimGenes); i++) {
    Gene *gene = Vector_getElementAt(prelimGenes, i);
    Transcript *transcript = Gene_getTranscriptAt(gene, 0);
    int giveUpFlag = 0; // Note: This is used a long way down this routine - if set results in skipping to next gene which is why I put it here
   
    // hack taking out weeny models

    if (Transcript_getLength(transcript) < 300) {
      continue;
    }

    Vector *models = Vector_new();

    int singleExon = 0;
    // first run on the fwd strand then on the reverse

//  STRAND: 
    int strand;
    for (strand = -1; strand<=1 && !giveUpFlag; strand+=2) {
      if ( RefineSolexaGenes_getRecursiveLimit(rsg) > 10000 ) {
        // set recursion to 10000 in case it was raised for a tricky gene
        RefineSolexaGenes_setRecursiveLimit(rsg, 10000);
        fprintf(stderr, "Warning: lowering recursive limit after complex gene\n"); 
      }

      fprintf(stderr,"Running on strand %d\n", strand);

// Doesn't seem to be used      StringHash *intronCount = StringHash_new(STRINGHASH_SMALL);
      Vector *exonIntron = Vector_new(); // A Vector of Vectors. Each Vector lists the possible introns for a particular exon
      
      StringHash *intronHash = StringHash_new(STRINGHASH_SMALL); // A hash containing all the introns, keyed on intron HitSeqName.
// Doesn't seem to be used      Vector *exonPrevIntron = Vector_new();
      StringHash *intronExon = StringHash_new(STRINGHASH_SMALL); // A Hash of Vectors. keyed on intron HitSeqName. Each element is a list of exons

      int mostRealIntrons = 0;
      double highestScore = 0;

      fprintf(stderr, "%s : %ld %ld:\n", Gene_getStableId(gene), Gene_getStart(gene), Gene_getEnd(gene));

      Vector *exons = RefineSolexaGenes_mergeExons(rsg, gene, strand);
      Vector_sort(exons, SeqFeature_startCompFunc);

/*
#      foreach my $exon ( @exons ) {
#        print  "EXTRAEXON: " . 
#          $exon->seq_region_name ." " .
#            ($exon->start +20) ." " .
#              ($exon->end -20)." " .
#                ( $exon->end - $exon->start -40)  ."\n"
#                  if $exon->{"_extra"} ;
#      }
*/
   
      int exonCount = Vector_getNumElement(exons);
      Vector *fakeIntrons = Vector_new();
      StringHash *knownExons = StringHash_new(STRINGHASH_SMALL);

      long offset = 0;

//    EXON:   
      int j;
      for (j=0; j<=exonCount; j++) {
        Exon *origExon = Vector_getElementAt(exons, j);
        Exon *exon = ExonUtils_cloneExon(origExon);
 
        int retainedIntron;
        int leftIntrons = 0;
        int rightIntrons = 0;

// Hack hack hack hack
// Don't seem to be used so commented out
//        $exon->{'left_mask'} = 0;
//        $exon->{'right_mask'} = $exon->length;


        printf("%d : %ld %ld\n",j, Exon_getStart(exon), Exon_getEnd(exon));
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

//      INTRON: 
        int k;
        for (k=0; k<Vector_getNumElement(introns); k++) {
          DNAAlignFeature *intron = Vector_getElementAt(introns, k);

          if (DNAAlignFeature_getStrand(intron) != strand) {
            continue;
          } 
          if (DNAAlignFeature_getLength(intron) > RefineSolexaGenes_getMinIntronSize(rsg)) {
            continue;
          } 
          if (DNAAlignFeature_getLength(intron) <= RefineSolexaGenes_getMaxIntronSize(rsg)) {
            continue;
          } 

          // discard introns that splice over our exon
          if (DNAAlignFeature_getStart(intron) < Exon_getStart(exon) && 
              DNAAlignFeature_getEnd(intron) > Exon_getEnd(exon) ) {
            intronOverlap++;
            continue;
          }

          // check to see if this exon contains a retained intron
          if (DNAAlignFeature_getStart(intron) > Exon_getStart(exon) && 
              DNAAlignFeature_getEnd(intron) < Exon_getEnd(exon) ) {
            retainedIntron = 1;
// Hack hack hack

            Exon_addFlag(exon, RSGEXON_RETAINED);

            Vector_addElement(retainedIntrons, intron);
          } else {
            // we need to know how many consensus introns we have to the 
            // left and right in order to determine whether to put in 
            // a non consensus intron
            if (DNAAlignFeature_getEnd(intron) <= Exon_getEnd(exon)) {
              if (strstr(DNAAlignFeature_getHitSeqName(intron), "non canonical")) {
                if (Intron_getScore(intron) > 1) {
                  Vector_addElement(leftNonConsIntrons, intron);
                }
              } else {
                Vector_addElement(leftConsIntrons, intron);
              }
            }
            if (DNAAlignFeature_getStart(intron) >= Exon_getStart(exon)) {
              if (strstr(DNAAlignFeature_getHitSeqName(intron), "non canonical")) {
                if (Intron_getScore(intron) > 1) {
                  Vector_addElement(rightNonConsIntrons, intron);
                }
              } else {
                Vector_addElement(rightConsIntrons, intron);
              }
            }
          }
        }
        
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

/* Does nothing
        if ( scalar(@left_c_introns)  == 0 && scalar(@left_nc_introns)  > 0) {
         # print STDERR "using " . scalar(@left_nc_introns) . " NC left \n";
        } 
        if ( scalar(@right_c_introns)  == 0 && scalar(@right_nc_introns)  > 0 ) {
         # print STDERR "using " . scalar(@right_nc_introns) . " NC right \n";
        }
*/
        
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
        Vector *exIntj = Vector_new();
        Vector_setElementAt(exonIntron, j, exIntj);

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

          StringHash_add(intronHash, DNAAlignFeature_getHitSeqName(intron), intron);

          // only use each intron twice once at the end and once at the start of
          // an exon
          // exon_intron links exons to the intron on their right ignoring strand
          if (DNAAlignFeature_getEnd(intron) > Exon_getEnd(exon)) {
            Vector_addElement(exIntj, intron);
          }

          // intron exon links introns to exons on their right ignoring strand
          if (DNAAlignFeature_getStart(intron) < Exon_getStart(exon)) {
            if (!StringHash_contains(intronExon, DNAAlignFeature_getHitSeqName(intron))) {
              StringHash_add(intronExon, DNAAlignFeature_getHitSeqName(intron), Vector_new());
            }
            Vector *ieVec = StringHash_getValue(intronExon, DNAAlignFeature_getHitSeqName(intron));
            Vector_addElement(ieVec, long_new(j));
            // exon_prev_intron links exons to introns on their left ignoring strand
// Doesn't seem to be used push @{ $exon_prev_intron[$j]}  , $intron ;
          }
        }
        if (Vector_getNumElement(retainedIntrons) > 0) {
          //#print STDERR "Dealing with " . scalar( @retained_introns ) . " retained introns \n";
          Vector *newExons = Vector_new();
          Vector_addElement(newExons, exon);

          // sort first by start then by end where start is the same
          Vector_sort(retainedIntrons, SeqFeature_startEndCompFunc);
/* This is just sorting by end after doing the start sort, so do that in one in the sort func
          int m;
// Think need the -1
          for (m=0; m < Vector_getNumElement(retainedIntrons)-1; m++) {
            DNAAlignFeature *retainedIntron = Vector_getElementAt(retainedIntrons, m);
            DNAAlignFeature *retainedIntronP1 = Vector_getElementAt(retainedIntrons, m+1);

            if (DNAAlignFeature_getStart(retainedIntron) == DNAAlignFeature_getStart(retainedIntronP1) &&
                DNAAlignFeature_getEnd(retainedIntron) > DNAAlignFeature_getEnd(retainedIntronP1)) {
              // reverse the order
              my $temp =  $retained_introns[m];
              $retained_introns[m] = $retained_introns[m+1];
              $retained_introns[m+1] = $temp;
            }
          }
*/

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
                //#print STDERR "RETAINED INTRON PENALTY for " . $intron->display_id ." before " . $intron->score . " ";
                rejectScore = DNAAlignFeature_getScore(intron) - RefineSolexaGenes_getRetainedIntronPenalty(rsg);
                // intron penalty is doubled for nc introns 
                if (strstr(DNAAlignFeature_getHitSeqName(intron), "non canonical")) {
                  rejectScore = rejectScore - RefineSolexaGenes_getRetainedIntronPenalty(rsg);
                }
                //#print STDERR " after " . $reject_score ."\n";
                if (rejectScore < 1 ) {
                  // treat as single exon
                  if (Vector_getNumElement(exons) == 1 ) {
                    // at least on this strand this model looks like a single exon
                    singleExon += 1;
                  }
                } else {
                  //#print STDERR  "Exon " . $exon->start ."\t". $exon->end . " has retained intron:\n     " . $intron->start ."\t" .  $intron->end ." "."\n";
                  // dont have circular references to exons or the paths
                  // will be infinite so clone this exon instead
                  // I guess we also want to keep the original exon too?
                  Exon *newExon1 = ExonUtils_cloneExon( exon );
                  Exon *newExon2 = ExonUtils_cloneExon( exon );
    
                  // chop it up a bit so it no longer overlaps the other introns
                  fprintf(stderr, "TRIMMING EXON\n");
                  int length = DNAAlignFeature_getEnd(intron) - DNAAlignFeature_getStart(intron);
                  Exon_setEnd  (newExon1, DNAAlignFeature_getStart(intron) + ( length / 2 ) - 2);
                  Exon_setStart(newExon2, DNAAlignFeature_getEnd(intron) - ( length / 2 ) + 2);
                  
                  char keKey[1024];
                  sprintf(keKey,"%ld-%ld-%d",Exon_getStart(newExon1), Exon_getEnd(newExon1), Exon_getStrand(newExon1));
                  if (!StringHash_contains(knownExons, keKey)) {
                    Vector_addElement(newExons, newExon1);
                    StringHash_add(knownExons, keKey, &trueVal);
                  }
  
                  sprintf(keKey,"%ld-%ld-%d",Exon_getStart(newExon2), Exon_getEnd(newExon2), Exon_getStrand(newExon2));
                  if (!StringHash_contains(knownExons, keKey)) {
                    Vector_addElement(newExons, newExon2);
                    StringHash_add(knownExons, keKey, &trueVal);
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
            Vector_removeElementAt(exons, i);
            for (k=0; k<Vector_getNumElement(newExons); k++) {
              Vector_insertElementAt(exons, i, Vector_getElementAt(newExons, k));
            }
/*  Doesn't seem to do anything
            for ( my $i = 0 ; $i<= $#exons ; $i++ ) {
              my $e = $exons[$i];
            }
*/
            printf("ADDED %d new exons\n", Vector_getNumElement(newExons));
            exonCount += Vector_getNumElement(newExons) -1; // was $#new_exons;
            // make sure they are all stil sorted
            Vector_sort(exons, SeqFeature_startCompFunc);
          }
        }
      }
      
      // replaced with if - next unless @exon_intron;

      if (Vector_getNumElement(exonIntron) > 0) {

        // Loop around the path generation, 
        // if there are too many paths to process return undef
        // then re-run the path processing but with increasing strictness
        // where strictness = elimianating alternate low scoring introns
        StringHash *paths;
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
          fprintf(stderr, "STRAND %d BEFORE COLLAPSING PATHS = %d\n", strand, StringHash_getNumValues(paths));
    
          char **pathKeys = StringHash_getKeys(paths);

          int nPath = StringHash_getNumValues(paths);

          // Perl sorted them 
          qsort(pathKeys, nPath, sizeof(char *), StrUtil_stringCompFunc);

          int j;
          // lets collapse redundant paths
          for (j=0; j<nPath; j++) {
            char *path = pathKeys[j];

            char **array;
            int nArray;

            //  print "PATHS $path\n";
            StrUtil_tokenizeByDelim(&array, &nArray, path, ".");

            char start[100000];  start[0] = '\0';
            char end[100000];    end[0] = '\0';
            char middle[100000]; middle[0] = '\0';
            

            int k;
            for (k=0; k<nArray; k++)  {
              //$start .= $array[k] ."."  unless k < 2;
              if (k>=2) sprintf(start,"%s%s.", start, array[k]);

              //$middle .= $array[k] ."." unless k < 2  or k >= $#array-1 ;
              if (k>=2 && k<nArray-2) sprintf(middle,"%s%s.",middle, array[k]);

              //$end .= $array[k] . "." unless k >= $#array-1;
              if (k<nArray-2) sprintf(end,"%s%s.", end, array[k]);
             
            }

            // remove redundancy from the array
            //delete $paths->{$start} if $start && $paths->{$start};
            if (start[0] != '\0' && StringHash_contains(paths, start)) {
              StringHash_remove(paths, start, NULL);
            }

            // delete $paths->{$end} if $end && $paths->{$end};
            if (end[0] != '\0' && StringHash_contains(paths, end)) {
              StringHash_remove(paths, end, NULL);
            }

            //delete $paths->{$middle} if $middle && $paths->{$middle};
            if (middle[0] != '\0' && StringHash_contains(paths, middle)) {
              StringHash_remove(paths, middle, NULL);
            }
          }
          

          for (j=0; j<nPath; j++) {
            free(pathKeys[j]);
          }
          free(pathKeys);
          
          fprintf(stderr, "AFTER COLLAPSING PATHS  = %d\n", StringHash_getNumValues(paths));
         
          Vector *strandModels = RefineSolexaGenes_makeModels(rsg, paths, strand, exons, gene, intronHash);
  
          Vector_append(models, strandModels);
          Vector_free(strandModels);
    
          fprintf(stderr, "Now have %d models\n", Vector_getNumElement(models));
        }
      }
    }

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
          clusteredModels = RefineSolexaGenes_reclusterModels(rsg, newClusters, &newClusters);
  //        print "Now have " .  scalar(@{$new_clusters}) ." new clusters after reclustering\n";
        }
      }
  
      if (clusteredModels != NULL && Vector_getNumElement(clusteredModels)) {
        Vector_append(cleanClusters, clusteredModels);
      }
      
      // filter to identify 'best', 'other' and 'bad' models
      RefineSolexaGenes_filterModels(rsg, cleanClusters);
  
      // process single exon models
      // if it has no introns on either strand
      if (RefineSolexaGenes_getSingleExonModelType(rsg) && singleExon == 2) {
        Vector *merged = RefineSolexaGenes_mergeExons(rsg, gene, 1);
        Exon *exon = Vector_getElementAt(merged, 0);
  
        Transcript *singleExonModel = NULL;
        //print STDERR " Single exon = $single_exon\n";
    
        if (Exon_getLength(exon) + 40 >= RefineSolexaGenes_getMinSingleExonLength(rsg)) {
          //print STDERR "Passed length filter " . $exon->length ."\n";
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
          Transcript *fwdTran = TranscriptUtils_computeTranslation(fwdT);
          
          long fwdTLen = 0;
          if (Transcript_getTranslation(fwdTran)) {
            fwdTLen = Translation_getGenomicEnd(Transcript_getTranslation(fwdTran)) - Translation_getGenomicStart(Transcript_getTranslation(fwdTran));
          }

          //print STDERR "FWD t length $fwd_t_len\n";
          Transcript *revT = Transcript_new();
          Transcript_addExon(revT, revExon, 0);
          //my $rev_t =  new Bio::EnsEMBL::Transcript(-EXONS => [$rev_exon]);
          
         // my $rev_tran = compute_translation(clone_Transcript($rev_t));
          Transcript *revTran = TranscriptUtils_computeTranslation(revT);
          
          long revTLen = 0;
          if (Transcript_getTranslation(revTran)) {
            revTLen = Translation_getGenomicEnd(Transcript_getTranslation(revTran)) - Translation_getGenomicStart(Transcript_getTranslation(revTran));
          }
    
          //print STDERR "REV t length $rev_t_len\n";
          if ( Transcript_getTranslation(fwdTran) &&  
               ( fwdTLen / Transcript_getLength(fwdTran)* 100 >= RefineSolexaGenes_getMinSingleExonCDSLength(rsg) &&
               fwdTLen >  revTLen )) {
            // keep this one
            singleExonModel =  fwdTran;
          }
    
          if ( Transcript_getTranslation(revTran) &&  
               ( revTLen / Transcript_getLength(revTran)* 100 >= RefineSolexaGenes_getMinSingleExonCDSLength(rsg) &&
               revTLen >  fwdTLen )) {
            // keep this one
            singleExonModel =  revTran;
          }
    
          if (singleExonModel != NULL) {
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
        }
      }
    }
// NIY: Any freeing/tidying for a gene's worth of processing should go here
  }
}

void RefineSolexaGenes_addToOutput(RefineSolexaGenes *rsg, Gene *gene) {
  if (rsg->output == NULL) {
    rsg->output = Vector_new();
  }
  Vector_addElement(rsg->output, gene);
}

Vector *RefineSolexaGenes_getOutput(RefineSolexaGenes *rsg) {
  return rsg->output;
}

// getAnalysis Should probably be in RunnableDB - put here for now
Analysis *RefineSolexaGenes_getAnalysis(RefineSolexaGenes *rsg) {
  return rsg->analysis;
}

Analysis *RefineSolexaGenes_createAnalysisObject(RefineSolexaGenes *rsg, char *logicName) {
  DBAdaptor *outDb = RefineSolexaGenes_getDbAdaptor(rsg, RefineSolexaGenes_getOutputDb(rsg));

  AnalysisAdaptor *aa = DBAdaptor_getAnalysisAdaptor(outDb);
  Analysis *analysis = AnalysisAdaptor_fetchByLogicName(aa, logicName);

  if (analysis != NULL) {
     return analysis ;
  }
  // need to create analysis object first
  analysis = Analysis_new();
  Analysis_setLogicName(analysis, logicName);

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
        if (Gene_getStrand(gene) == strand) {
        //#print "GENE STRAND " . $gene->strand ." vs $strand \n";
          Vector_addElement(strandedGenes, gene);
        }
      }

      if (Vector_getNumElement(strandedGenes) == 0) {
        Vector_free(strandedGenes);
        continue;
      }

      ModelCluster *strandedCluster = RefineSolexaGenes_recalculateCluster(rsg, strandedGenes);

      Gene *best = NULL;
      Vector *genes;
      Vector *otherGenes;

      for (j=0; j<Vector_getNumElement(strandedCluster->finalModels); j++) {
        Gene *gene = Vector_getElementAt(strandedCluster->finalModels, j);

        if (Gene_getStrand(gene) == strand) {
          if (!strcmp(Gene_getBiotype(gene), RefineSolexaGenes_getBestScoreType(rsg))) {
            //$best =  $gene if $gene->biotype eq $self->BEST_SCORE;
            best = gene;
            Vector_addElement(genes, gene);
          }
        }
      }

      if (best == NULL) {
        fprintf(stderr,"No best model found\n");
        exit(1);
      }

      // now recluster 
    //OTHERGENE: 
      for (j=0; j<Vector_getNumElement(strandedCluster->finalModels); j++) {
        Gene *gene = Vector_getElementAt(strandedCluster->finalModels, j);

        if (!strcmp(Gene_getBiotype(gene), RefineSolexaGenes_getBestScoreType(rsg))) {
          continue;
        }

        Transcript *bestTrans = Gene_getTranscriptAt(best, 0);
        Vector *bestExons = Vector_copy(Transcript_getAllExons(bestTrans));
        Vector_sort(bestExons, SeqFeature_startCompFunc);

        Transcript *geneTrans = Gene_getTranscriptAt(gene, 0);
        Vector *otherExons = Vector_copy(Transcript_getAllExons(geneTrans));
        Vector_sort(otherExons, SeqFeature_startCompFunc);

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
                
             //   print "Overlap " . $gene->start . " " , $gene->end . " " . $gene->strand ."\n";
  // Big jump here - think about how to do
                //next OTHERGENE;
                doneBest = 1;
                break;
              }
            }
          }
        }
        // other model has no exon overlap with best model it needs to be in a new cluster
        if (!doneBest) {
          Vector_addElement(otherGenes, gene);
        }
        //#print "No overlap " . $gene->start . " " , $gene->end . " " . $gene->strand ."\n";
      }

      // now we need to fix the clusters
      if (Vector_getNumElement(otherGenes) > 0) {
        Vector_addElement(finalClusters, RefineSolexaGenes_recalculateCluster(rsg, genes));
        Vector_addElement(newClusters,   RefineSolexaGenes_recalculateCluster(rsg, otherGenes));
// NIY: Free strandedCluster???
      } else {
        // keep it as it was
        Vector_addElement(finalClusters, strandedCluster);
      }
// NIY: Free stuff???

    }
  }

  // print "CLUSTERS " . scalar(@final_clusters) ." " .  scalar(@new_clusters) ."\n";
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

    Transcript *t = Gene_getTranscriptAt(g, i);
    if (Transcript_getScore(t) > score) {
      score = Transcript_getScore(t);
    }
  }

  cluster->start  = start;
  cluster->end    = end;
  cluster->strand = strand;

  for (i=0; i<Vector_getNumElement(genes); i++) {
    Gene *g = Vector_getElementAt(genes, i);
    Transcript *t = Gene_getTranscriptAt(g, i);

    if (Transcript_getScore(t) == score ) {
      Gene_setBiotype(g, RefineSolexaGenes_getBestScoreType(rsg));
      // print "BEST SCORE " . $g->start ." " . $g->end . " " . $g->strand . "\n";
    } else {
      Gene_setBiotype(g, RefineSolexaGenes_getOtherIsoformsType(rsg));
    }
  }

  cluster->finalModels = genes;

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

  Vector *models = Vector_new();

  int i;
  for (i=0; i<Vector_getNumElement(clusters); i++) {
    ModelCluster *cluster = Vector_getElementAt(clusters, i);

    if (cluster->finalModels == NULL) {
      fprintf(stderr,"Had a cluster with no final models\n");
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

  // overlaps
  for (i=0; i<Vector_getNumElement(fwd); i++) {
    ModelCluster *fc = Vector_getElementAt(fwd, i);
    int j;
    for (j=0; j<Vector_getNumElement(rev); j++) {
      ModelCluster *rc = Vector_getElementAt(fwd, i);

      // one is within the other  or they are the same
      // they proably need to be rejected on the basis of coding overlap
      if (( fc->start >= rc->start && fc->end <= rc->end) || 
          ( rc->start >= fc->start && rc->end <= fc->end))  {
        
        // do they have coding overlap?
        Vector *fgs = fc->finalModels;
        Vector *rgs = rc->finalModels;

        // do they have coding overlap?
//      FG: 
        int k;
        for (k=0; k<Vector_getNumElement(fgs); k++) {
          Gene *fg = Vector_getElementAt(fgs, k);
          Transcript *ft = Gene_getTranscriptAt(fg, 0);
          int done = 0;
          
          if (!Transcript_getTranslation(ft)) {
            continue;
          }

          if (Translation_getLength(Transcript_getTranslation(ft)) <= 100) {
            Gene_setBiotype(fg, "bad");
            // Do an else instead next FG;
          } else {

            int m;
// Maybe leak here - translateableExons can be allocated
            Vector *ftTranslateableExons = Transcript_getAllTranslateableExons(ft);
            for (m=0; m<Vector_getNumElement(ftTranslateableExons) && !done; m++) {
              Exon *fe = Vector_getElementAt(ftTranslateableExons,  m);

//            RG: 
              int n;
              for (n=0; n<Vector_getNumElement(rgs) && !done; n++) {
                Gene *rg = Vector_getElementAt(rgs, n);
                Transcript *rt = Gene_getTranscriptAt(rg, 0);

// Next what?????
                if (!Transcript_getTranslation(rt)) {
                  continue;
                }
                if (Translation_getLength(Transcript_getTranslation(rt)) <=  100 ) {
                  Gene_setBiotype(rg, "bad");
                  // Do an else instead next RG;
                } else {
                  int p;
                  Vector *rtTranslateableExons = Transcript_getAllTranslateableExons(rt);
                  for (p=0; p<Vector_getNumElement(rtTranslateableExons) && !done; p++) {
                    Exon *re = Vector_getElementAt(rtTranslateableExons, p);

                    if ( Exon_getStart(fe) <= Exon_getEnd(re) && 
                         Exon_getEnd(fe)  >=  Exon_getStart(re)) {
                      // coding overlap        

// Hack hack hack
                      if (Transcript_getScore(ft) < Transcript_getScore(rt)) {
                        //# get rid of / label the reverse genes 
                        Gene_setBiotype(fg, "bad");
                      } else {
                        Gene_setBiotype(rg, "bad");
                      }
// Big jump - think about how to do
                      //next FG;
                      done = 1;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

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

      //print "TRANSCRIPT " . $gene->get_all_Transcripts->[0]->{'_depth'} .           " Exon use $exon_use Biotype " . $gene->biotype ."\n";
      int es = 0;
      int ee = 0;
      char pattern[65500];
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
        //        print "CALLING it bad\n";
      }

      if (StringHash_contains(exonPattern, pattern)) {
        // seen it before - or something very much like it
        Gene_setBiotype(gene, "duplicate");
        //print "CALLING it bad\n";
      } else {
        StringHash_add(exonPattern, pattern, &trueVal);
      }
    } 

    // promote "bad" models that have a cds as long as the best cds to 
    // alt isoforms
    Vector *finalModels = cluster->finalModels;
    int bestCds = 0;
    int g;
    for (g=0; g<Vector_getNumElement(finalModels); g++) {
      Gene *gene = Vector_getElementAt(finalModels, g);

      Transcript *transcript = Gene_getTranscriptAt(gene, 0);
      
      fprintf(stdout, "%d - %f tran length %d\n", 
              g, Transcript_getScore(transcript), Transcript_getcDNACodingEnd(transcript) - Transcript_getcDNACodingStart(transcript));

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
          }
        }
      } else {
        // if NOT duplicate 
        if (strcmp(Gene_getBiotype(gene), "duplicate")) {
          if (!strcmp(Gene_getBiotype(gene), RefineSolexaGenes_getBestScoreType(rsg))) {
            // trim the UTR
            fprintf(stderr,"UTR trimming not implemented yet\n");
// NIY            $self->prune_UTR($gene);
            RefineSolexaGenes_addToOutput(rsg, gene);
          } else {
// Note here checking for whether other isoforms type is null, other places assume its not null
            if (RefineSolexaGenes_getOtherNum(rsg) && RefineSolexaGenes_getOtherIsoformsType(rsg) != NULL && count <= RefineSolexaGenes_getOtherNum(rsg)) {
              // trim the UTR
              fprintf(stderr,"UTR trimming not implemented yet\n");
// NIY              $self->prune_UTR($gene);
              RefineSolexaGenes_addToOutput(rsg, gene);
            }
          }
        }
      }
      if (!strcmp(Gene_getBiotype(gene), RefineSolexaGenes_getOtherIsoformsType(rsg)) || 
          !strcmp(Gene_getBiotype(gene), RefineSolexaGenes_getBestScoreType(rsg))) {
        count++;
      }
    }
  }
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
  Vector *genes = Vector_new();

  char **pathKeys = StringHash_getKeys(paths);

  int i;
  for (i=0; i<StringHash_getNumValues(paths); i++) {
    char *path = pathKeys[i];

    char exonUse[8000]; exonUse[0] = '\0';

    Model *model;
    if ((model = (Model *)calloc(1,sizeof(Model))) == NULL) {
      fprintf(stderr,"Failed allocating model\n");
      exit(1);
    } 

    model->features = Vector_new();
    double exonScore = 0;
    double intronScore = 0;

    char **featureStrings;
    int nFeatureString;

    StrUtil_tokenizeByDelim(&featureStrings, &nFeatureString, path, ".");

    int j;
    for (j=0; j<nFeatureString; j++)  {
      char *featureString = featureStrings[j];
      if (strstr(featureString, "canonical")) { // Intron
        DNAAlignFeature *intron = StringHash_getValue(intronHash, featureString);
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

    double totalScore = ((int)exonScore)/100 + intronScore;
    // last elements are the strand and score
    
    StrUtil_copyString(&model->exonUse, exonUse, 0);
    model->totalScore = totalScore;
    Vector_addElement(models, model);

    // NIY: Free stuff
  }

  printf("Starting model_cluster\n");
  // now lets cluster the models so that they are non overlapping
  // and return the clusters arranged by score

  Vector *modelClusters = RefineSolexaGenes_makeModelClusters(rsg, models, strand);
  
  printf("Starting gene cycle\n");
  // Now we cycle through all the models and turn them into genes
  // we start with the highest scoring modes and work backwards
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
      double intronScore = 0;
      double exonScore = 0;
      int nonConIntrons = 0;

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
            exit(1);
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
          if (strstr(DNAAlignFeature_getHitSeqName(intron), "non canonical")) {
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

          if (strstr(DNAAlignFeature_getHitSeqName(intron), "non canonical")) {
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
        //next MODEL;
// NIY: Free stuff??
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
// NIY: Free stuff??
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
        continue;
      }

      // make it into a gene
      Vector *modifiedExons = Vector_new();
      for (k=0; k<Vector_getNumElement(newExons); k++) {
        SeqFeature *sf = Vector_getElementAt(newExons, k);
        if (sf->objectType != CLASS_DNADNAALIGNFEATURE) {
// Not sure if I really need to make these clones
          Vector_addElement(modifiedExons, ExonUtils_cloneExon((Exon *)sf));
        }
      }
 
      if ( strand == 1 ) {
        Vector_sort(modifiedExons, SeqFeature_startCompFunc);
      } else {
        Vector_sort(modifiedExons, SeqFeature_reverseStartCompFunc);
      }

      // make it into a gene
//      my $t =  new Bio::EnsEMBL::Transcript(-EXONS => \@modified_exons);
      Transcript *t = Transcript_new();
      for (k=0; k<Vector_getNumElement(modifiedExons); k++) {
        Transcript_addExon(t, (Exon *)Vector_getElementAt(modifiedExons, k), k);
      }
 
      for (k=0; k<Vector_getNumElement(ises); k++) {
        IntronSupportingEvidence *ise = Vector_getElementAt(ises, k);
        Transcript_addIntronSupportingEvidence(t, ise);
      }
// NIY: Free ises vector??

//# check for dna
//#      my $check = $t->seq->seq ;
//#      my $Ns =  $check =~  s/N//g;
//#      if( length($t->seq->seq) == $Ns ){
//#        $self->throw("There does not appear to be ay DNA in the database, transcript seq is all N's\n");
//#      }
      // add a translation 
      // NOTE: Extra cloning my $tran = compute_translation(clone_Transcript($t));
      Transcript *tran = TranscriptUtils_computeTranslation(t);

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
            DNAAlignFeature *daf = Vector_getElementAt(support, 0);

            codingBonus += DNAAlignFeature_getScore(daf);
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
      Transcript_setScore(tran, (int)(( intronScore + exonScore ) / 10 ) + codingBonus);

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
    }

    // re-sort the transcripts to take account of the revised scores
    Vector_sort(trans, SeqFeature_reverseScoreCompFunc);
    //@trans = sort { $b->{'_score'} <=> $a->{'_score'} } @trans;

// Doesn't seem to be used (only set)    Transcript *best;
    for (j=0; j<Vector_getNumElement(trans); j++) {
      Transcript *tran = Vector_getElementAt(trans, j);

      TranscriptExtraData *ted = Transcript_getExtraData(tran);

      Gene *newGene = TranscriptUtils_convertToGene(tran, Gene_getAnalysis(gene), NULL);

      version++;
      Gene_setBiotype(newGene, RefineSolexaGenes_getOtherIsoformsType(rsg));

      if (version == 1) {
        Gene_setBiotype(newGene, RefineSolexaGenes_getBestScoreType(rsg));
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
  }

  printf("Done gene cycle\n");
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
  Gene_setBiotype   (g, Transcript_getBiotype(t));
  Gene_addTranscript(g, Transcript_getBiotype(t));
  Gene_setAnalysis  (g, analysis);
// Odd!!!
  Gene_setDbID      (g, Transcript_getDbID(t));
  
  
/* Wierd checks
    throw("there are no transcripts for this gene\n") if scalar(@{$g->get_all_Transcripts}) == 0 ;
    for my $tr ( @{$g->get_all_Transcripts} ) {
      throw("there are no exons  for this transcript \n") if scalar(@{$tr->get_all_Exons}) == 0 ;
    }
    throw("Problems with ".Gene_info($gene)." undef coords") if(!$gene->start || !$gene->end);
*/

  return g;
}


/* NIY
Gene *RefineSolexaGenes_pruneUTR(RefineSolexaGenes *rsg, Gene *gene) {

  if (!RefineSolexaGenes_trimUTR(rsg)) {
    return gene;
  }

  my $transcript = $gene->get_all_Transcripts->[0];
  unless ( $transcript->translateable_seq ) {
    return $gene;
  }
 
  // fetch introns 
  my $introns = $transcript->{'_introns'};

  // otherwise trim the UTR according to the values set out in the config
  my %intron_hash;

  foreach my $intron ( @{$introns} ) {
    my $key = $intron->start  .":". $intron->end .":". $intron->strand;
    $intron_hash{$key} = $intron;
  }
  
  my @new_fivep;
  my @new_threep;
  my @new_exons;
  my @features;
  my @exons = sort {$a->start <=> $b->start }  @{$transcript->get_all_Exons};
  
  // put everything into the features array
  push @features, @exons;
  for ( my $i =0 ; $i < $#exons  ; $i++ ) {
    my $key = ($exons[$i]->end) .":". ($exons[$i+1]->start ) . ":" . $exons[$i]->strand;    
    if ( my $intron = $intron_hash{$key}  ) {
      push @features, $intron;
    }
  }
  @features = sort { $a->start <=> $b->start } @features;
  // so now we should have an array of alternating introns and exons
  print "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Trimming UTR
Transcript " .  $transcript->seq_region_name ." " . 
    $transcript->start ." " . 
      $transcript->end ." " . 
        $transcript->strand ." " . 
          scalar(@{$transcript->get_all_Exons}) ."

";
  throw("Something is wrong we are missing " . scalar(@{$introns}) ." introns " . scalar(@exons) . "  exons " . scalar(@features) . " exons and introns\n")
    unless scalar(@features) == (scalar(@exons) * 2) -1 ;
  my $average_intron = 0;
  my $intron_count = 0;
  // leave single exon genes alone for now
  if ( scalar(@features) == 1 or scalar(@{$transcript->get_all_translateable_Exons}) == 1 )  {
    // lets strip the UTR
    my $trimmed_tran =  $self->modify_transcript($transcript,$transcript->get_all_translateable_Exons);
    // The naughty bit!
    $gene->{_transcript_array} = [];
    $gene->add_Transcript($trimmed_tran);
    return   $gene;
  }
  // first calculate the average
  foreach my $f ( @features ) {
    if ( $f->isa("Bio::EnsEMBL::DnaDnaAlignFeature")) {
      $average_intron += $f->score;
      $intron_count++;
    }
  }
  $average_intron /= $intron_count;

  foreach my $f ( @features ) {
   //# print $f->start . "\t" . $f->end ."\t$f\t";
    if ( $f->isa("Bio::EnsEMBL::DnaDnaAlignFeature")) {
    //#  print $f->score;
      if ( $average_intron && ($f->score /  $average_intron) * 100 <= $self->REJECT_INTRON_CUTOFF ) {
        print " Potentially bad ";
      }
    }
  //#  print "\n";
  }
  throw("Something is wrong we are missing introns " . scalar(@exons) . "  exons  and $intron_count introns\n")
    unless $intron_count == scalar(@exons) -1 ;
  print  "Average intron depth = $average_intron \n";
  
  
  my @fivep;
  my @threep;
  my $coding =0 ;
  // need to account for strand
  @features = sort { $b->start <=> $a->start } @features if $transcript->strand == -1;

  for ( my $i = 0 ; $i <=  $#features ; $i += 2  ) {
    my $e = $features[$i];
    throw("Got a DNA align feature where I should have an exon\n")
      unless $e->isa("Bio::EnsEMBL::Exon");
    if ( $e->coding_region_start($transcript) ) {
      // first coding exon
      for ( my $j = 0 ; $j <= $i ; $j++ ) {
        push @fivep,$features[$j];
      }
      last;
    }
  }
  for ( my $i = $#features ; $i >= 0 ;  $i -= 2  ) {
    my $e = $features[$i];
    throw("Got a DNA align feature where I should have an exon\n")
      unless $e->isa("Bio::EnsEMBL::Exon");
    if ( $e->coding_region_end($transcript) ) {
          # first coding exon
      for ( my $j = $i ; $j <= $#features ; $j++ ) {
        push @threep,$features[$j];
      }
      last;
    }
  }
  
  // want to start at last coding exon and work outwards so....
  @fivep = reverse @fivep;
  // now we should be good
  print "FIVE P \n";
  @new_exons = @{$transcript->get_all_translateable_Exons};
  my $fivep_cds = shift(@new_exons);
  my $threep_cds = pop(@new_exons);
  my $fiveplen;
  my $threeplen;
  my $fivepc = 0 ;
  my $threepc = 0 ;
  my $nmd;
  
  // FIVE PRIME RULES
  
 FIVEP: for ( my $i = 0 ; $i <= $#fivep ; $i++ ) {
    my $f =  $fivep[$i];
    print $f->start . "\t" . $f->end ."\t$f\t";
    if ( $i == 0 ) {
      throw("First feature is not an exon \n")
        unless $f->isa("Bio::EnsEMBL::Exon");
      // UTR starts in this exon - how long is it?
      my $cds_start = $f->coding_region_start($transcript);
      $cds_start = $f->coding_region_end($transcript)  if $transcript->strand == -1;
      throw("First coding exon has no CDS \n") unless $cds_start;
      print "CDS START $cds_start\t";
      $fiveplen = $cds_start - $f->start +1 if $transcript->strand == 1;
      $fiveplen = $f->end - $cds_start   +1 if $transcript->strand == -1;
      // is the coding exon too long
      if ( $fiveplen > $self->MAX_5PRIME_LENGTH ) {
        // replace it with the cds
        @new_fivep = ($fivep_cds);
        print " 5p too long $fiveplen \n";
        last FIVEP;
      }
      push @new_fivep,$f;
      $fivepc++;
    } else {
      if (  $f->isa("Bio::EnsEMBL::Exon") ) {
        $fivepc++;
        $fiveplen+= $f->end - $f->start +1;
        // does it make the UTR too long?
        if ( $fiveplen > $self->MAX_5PRIME_LENGTH ) {
          // dont add it
          print " 5p too long $fiveplen \n";
          last FIVEP;
        }
        // is it too many exons?
        if ( $fivepc > $self->MAX_5PRIME_EXONS ) {
          // dont add it
          print " too many 5p  $fivepc cut them all as we are not sure \n";
          @new_fivep = ($fivep_cds);
          last FIVEP;
        }
        push @new_fivep,$f;
      }
    }
    // Does the intron score well enough to include the exon
    // apply rules and add successful exons into the mix
    if ( $f->isa("Bio::EnsEMBL::DnaDnaAlignFeature")) {
      print $f->score;
      if ( $average_intron && ($f->score /  $average_intron) *100 <= $self->REJECT_INTRON_CUTOFF ) {
        print " Rejecting on score cutoff " . $f->score ." vs  $average_intron\n";
        // dont add any more 
        last FIVEP;
      }
    }
    print "\n";
  }
  
  // three P
  print "THREE P \n";
 THREEP:   for ( my $i = 0 ; $i <= $#threep ; $i++ ) {
    my $f = $threep[$i];
    print $f->start . "\t" . $f->end ."\t$f\t";
    if ( $i == 0 ) {
      throw("First feature is not an exon \n")
        unless $f->isa("Bio::EnsEMBL::Exon");
      // UTR starts in this exon - how long is it?
      my $cds_end = $f->coding_region_end($transcript);
      $cds_end = $f->coding_region_start($transcript)  if $transcript->strand == -1;
      throw("last coding exon has no CDS \n") unless $cds_end;
      print "CDS END $cds_end\t";
      $threeplen = $cds_end - $f->start +1 if $transcript->strand == -1;
      $threeplen = $f->end - $cds_end   +1 if $transcript->strand == 1;
      # is the coding exon too long
      if ( $threeplen > $self->MAX_3PRIME_LENGTH ) {
        // replace it with the cds
        @new_threep = ($threep_cds);
        print " 3p too long $threeplen \n";
        last THREEP;
      }
      push @new_threep,$f;
      $nmd = $threeplen ;
      $threepc++;
    } else {
      if (  $f->isa("Bio::EnsEMBL::Exon") ) {
        // does it break the NMD rule?
        if ( $nmd > 55 ) {
          print " splice is after $nmd bp from stop codon - rejected on NMD rule of maximum 55 bp \n";
          @new_threep = ($threep_cds);
          last THREEP;
        }
        $threepc++;
        $threeplen+= $f->end - $f->start +1;
        // does it make the UTR too long?
        if ( $threeplen > $self->MAX_3PRIME_LENGTH ) {
          // dont add it
          print " 3p too long $threeplen \n";
          last THREEP;
        }
        // is it too many exons?
        if ( $threepc > $self->MAX_3PRIME_EXONS ) {
          // dont add it
          print " too many 3p  $threepc cut them all as we are not sure \n";
          @new_threep = ($threep_cds);
          last THREEP;
        }
        push @new_threep,$f;
      }
    }
    // Does the intron score well enough to include the exon
    // apply rules and add successful exons into the mix
    if ( $f->isa("Bio::EnsEMBL::DnaDnaAlignFeature")) {
      print $f->score;
      if ( $average_intron && ($f->score /  $average_intron) * 100 <= $self->REJECT_INTRON_CUTOFF ) {
        print " Rejecting on score cutoff " . $f->score ." vs  $average_intron\n";
        // dont add any more 
        last THREEP;
      }
    }
    print "\n";
  }
  
  push @new_exons, @new_fivep;
  push @new_exons, @new_threep;
  print " New transript has " . scalar(@new_exons) , " exons\n";
  my @clones;
  foreach my $e ( @new_exons ) {
    throw("Not is not an exon " . $e->start ." " . $e->end . " $e\n")
      unless $e->isa("Bio::EnsEMBL::Exon");
    push @clones, clone_Exon($e);
  }
  @clones = sort { $a->start <=> $b->start } @clones;
  @clones =  reverse(@clones) if $transcript->strand == -1;
  my $trimmed_tran =  $self->modify_transcript($transcript,\@clones);
  // The naughty bit!
  $gene->{_transcript_array} = [];
  $gene->add_Transcript($trimmed_tran);
  return $gene;
}
*/


// NOTE: Exons MUST BE sorted in correct order for adding to transcript
Transcript *RefineSolexaGenes_modifyTranscript(RefineSolexaGenes *rsg, Transcript *tran, Vector *exons) {
  long cdsStart;
  long cdsEnd;

  if (Transcript_getStrand(tran) == 1) {
    cdsStart = Transcript_getCodingRegionStart(tran);
    cdsEnd   = Transcript_getCodingRegionEnd(tran);
  } else {
    cdsEnd   = Transcript_getCodingRegionStart(tran);
    cdsStart = Transcript_getCodingRegionEnd(tran);
  }
  
  printf("CDS START END %ld %ld\n", cdsStart,  cdsEnd);
  Translation *tln = Transcript_getTranslation(tran); 
  printf("PHASE %d %d\n", Translation_getStart(tln), Translation_getEnd(tln));

// t created here and again below - not sure if need both
  Transcript *t = Transcript_new();
  int i;
  for (i=0; i<Vector_getNumElement(exons); i++) {
    Exon *exon = Vector_getElementAt(exons, i);
    Transcript_addExon(t, exon, i);
  }

  Exon *se;
  Exon *ee;
  for (i=0; i<Transcript_getExonCount(t); i++) {
    Exon *e = Transcript_getExonAt(t, i);
    if (Exon_getStart(e) <= cdsEnd && Exon_getEnd(e) >= cdsEnd) {
      ee = e;
    }
    if (Exon_getStart(e) <= cdsStart && Exon_getEnd(e) >= cdsStart) {
      se = e;
    }
  }

  long ts;
  long te;
  if ( Transcript_getStrand(tran) == -1 ) {
    ts =  Exon_getEnd(se) - cdsStart + 1;
    te =  Exon_getEnd(ee) - cdsEnd   + 1;
  } else {
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
      }
    }
  }
  //#  my $start_phase = $se->phase;
  Translation *translation = Translation_new();
  Translation_setStartExon(translation, se); 
  Translation_setEndExon(translation, se); 
  Translation_setStart(translation, ts); 
  Translation_setEnd(translation, te); 

  printf("S-E %ld %ld\n",ts, te); //#START PHASE $start_phase\n";
  printf("GS %ld %ld\n", Translation_getGenomicStart(translation), Translation_getGenomicEnd(translation));
  Transcript_setTranslation(t, translation);
  //# calculate_exon_phases($t,$start_phase);

  char *tranTranslationSeq = Transcript_translate(tran);
  char *tTranslationSeq = Transcript_translate(t);
  if (strcmp(tranTranslationSeq, tTranslationSeq)) {
    fprintf(stderr, "Translations do not match: Before %s\nAfter  %s\n", Transcript_translate(tran), Transcript_translate(t));
    exit(1);
  }
// NIY: free translation strings??
  return t;
}

void RefineSolexaGenes_writeOutput(RefineSolexaGenes *rsg) {
  DBAdaptor *outdb = RefineSolexaGenes_getDbAdaptor(rsg, RefineSolexaGenes_getOutputDb(rsg));

  GeneAdaptor *geneAdaptor = DBAdaptor_getGeneAdaptor(outdb);
// NIY:
//  $outdb->dbc->disconnect_when_inactive(0);

  Vector *output = RefineSolexaGenes_getOutput(rsg);
  
  int fails = 0;
  int total = 0;
//  GENE: 
  int i;
  for (i=0; i<Vector_getNumElement(output); i++) {
    Gene *gene = Vector_getElementAt(output, i);
    Analysis *anal = RefineSolexaGenes_getAnalysis(rsg);

    Gene_setAnalysis(gene, anal);
    Gene_setSource(gene, Analysis_getLogicName(anal));

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
          continue; // next GENE
        }
        Gene_setBiotype(gene, RefineSolexaGenes_getSingleExonModelType(rsg));
      } else {
        // dont store it
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
  }

  DNAAlignFeatureAdaptor *intronAdaptor = DBAdaptor_getDNAAlignFeatureAdaptor(outdb);
  fails = 0;
  total = 0;
 
  Vector *intronFeatures = RefineSolexaGenes_getIntronFeatures(rsg);

  for (i=0; i<Vector_getNumElement(intronFeatures); i++) {
    DNAAlignFeature *intron = Vector_getElementAt(intronFeatures, i);

//# SMJS For now leave as exon end and exon start so that it edge matches in apollo
//#    $intron->start($intron->start+1);
//#    $intron->end($intron->end-1);
//    eval {
        DNAAlignFeatureAdaptor_store(intronAdaptor, intron);
//    };
/*
    if ($@){
      warning("Unable to store DnaAlignFeature!!\n$@");
      $fails++;
    }
*/
    total++;
  }
  if (fails > 0) {
    fprintf(stderr, "Not all introns could be written successfully (%d fails out of %d)\n", fails, total);
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

#define PT_ERROR 1
int RefineSolexaGenes_processTree(RefineSolexaGenes *rsg, StringHash *hashref, char *index, char *soFar, StringHash *paths) {
  // dont let it go on for ever eating memory
  if (limit > RefineSolexaGenes_getRecursiveLimit(rsg)){
    fprintf(stderr,"Too many recursive possibilities\n");
    return PT_ERROR;
  }

// Perl passes sofar by value so makes lots of strings with the different paths
// To replicate that, make copy of sofar here - must be allocated off heap, pointer to it is on stack, so should emulate the perl way
  char *newSoFar = StrUtil_copyString(&newSoFar, soFar, 0);
  char tmpStr[512];
  sprintf(tmpStr,"%s.",index);
  newSoFar = StrUtil_appendString(newSoFar, tmpStr);

  int nNode = 0;
  if (StringHash_contains(hashref, index)) {
    StringHash *subHash = StringHash_getValue(hashref, index);
  
    char **node = StringHash_getKeys(subHash);
    nNode = StringHash_getNumValues(subHash);

    int i;
    for (i=0; i<nNode; i++) {
      char *child = node[i];

      limit++;
      int result = RefineSolexaGenes_processTree(rsg, hashref, child, newSoFar, paths);
      if (result == PT_ERROR) {
        limit = 0;
        return PT_ERROR; 
      }
    }
  }
  if (nNode == 0) {
    printf("path key when nNode = 0 is %s\n", soFar);
    StringHash_add(paths, newSoFar, &trueVal);
  }
  return paths;
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
  StringHash *variants;
  int removed = 0;
  int i;
  int j;

  //# now lets make a hash of hashes holding which exons connect to which
  for (i=0; i<Vector_getNumElement(exons); i++) {
    Exon *exon = Vector_getElementAt(exons, i);
    Vector *exInti = Vector_getElementAt(exonIntron, i);
    if (Vector_getNumElement(exInti)) {
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
        Vector *intronGroup;
        for (j=0; j<Vector_getNumElement(exInti); j++) {
          DNAAlignFeature *intron = Vector_getElementAt(exInti, j);

          if (strstr(DNAAlignFeature_getHitSeqName(intron), "REMOVED")) {
            Vector_addElement(intronGroup, intron);
          }
        }

        //# now lets sort these groups by score
        Vector_sort(intronGroup, SeqFeature_reverseScoreCompFunc);

        //# now lets see what they look like
        fprintf(stderr,"EXON %d: %ld - %ld - %d\n",
                i, Exon_getStart(exon), Exon_getEnd(exon), Exon_getStrand(exon));

        for (j=0; j<Vector_getNumElement(intronGroup); j++) {
          DNAAlignFeature *intron = Vector_getElementAt(intronGroup, j);
          fprintf(stderr, "%d %s %f\n", i, DNAAlignFeature_getHitSeqName(intron), DNAAlignFeature_getScore(intron));
        }

        if (Vector_getNumElement(intronGroup) > 1) {
          // remove the lowest scoring one
          //my $intron = pop( @{$intron_groups{$group}} ) ;
          DNAAlignFeature *intron = Vector_getLastElement(intronGroup);
          Vector_removeElementAt(intronGroup, Vector_getNumElement(intronGroup)-1);

          fprintf(stderr,"Eliminating %s %f\n",DNAAlignFeature_getHitSeqName(intron), DNAAlignFeature_getScore(intron));
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

              //$variants->{$i}->{$intron->hseqname} = 1;
              if (!StringHash_contains(variants, exonIndStr)) {
                StringHash_add(variants, exonIndStr, StringHash_new(STRINGHASH_SMALL));
              }
              StringHash *varEHash = StringHash_getValue(variants, exonIndStr);
              StringHash_add(varEHash, hseqname, &trueVal);

              //$variants->{$intron->hseqname}->{$exon} = 1;
              if (!StringHash_contains(variants, hseqname)) {
                StringHash_add(variants, hseqname, StringHash_new(STRINGHASH_SMALL));
              }
              StringHash *varIHash = StringHash_getValue(variants, hseqname);
              StringHash_add(varIHash, exonIndStr, &trueVal);
            }
          }
        }
      }
    }
  }
  
// ??? why &! rather than && ! if ($strict &! $removed ) {
  if (strict && !removed) {
    Exon *firstExon = Vector_getElementAt(exons, 0);

    fprintf(stderr,"Warning: Cannot simplify this gene any more EXON 0: %ld - %ld - %d\n",
            Exon_getStart(firstExon), Exon_getEnd(firstExon), Exon_getStrand(firstExon));

    if (RefineSolexaGenes_getRecursiveLimit(rsg) < RefineSolexaGenes_getMaxRecursions(rsg)) {
      RefineSolexaGenes_setRecursiveLimit(rsg, RefineSolexaGenes_getRecursiveLimit(rsg) * 10);
      fprintf(stderr, "Upping recursive limit to %d to see if it helps\n", RefineSolexaGenes_getRecursiveLimit(rsg));
    } else {
      fprintf(stderr,"Giving up on EXON 0: %ld - %ld - %d\n",
              Exon_getStart(firstExon), Exon_getEnd(firstExon), Exon_getStrand(firstExon));
//!!!!!!!!!!! NIY What to return
      *giveUpFlag = 1;
      return NULL;
    }
  }
  
  // work out all the possible paths given the features we have
  int result;
  StringHash *paths = StringHash_new(STRINGHASH_LARGE);
  for (i=0; i<Vector_getNumElement(exons) && result != PT_ERROR; i++) {
    limit = 0;
// Instead of undef for sofar I pass in ""
    char tmpStr[50];
    sprintf(tmpStr,"%d",i);
    result = RefineSolexaGenes_processTree(rsg, variants, tmpStr, "", paths);

    if (result == PT_ERROR) {
      fprintf(stderr, "Could not process cluster %d trying again with simpler cluster\n", i);
      return NULL;
    }
  }

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

  // sort them by the start of the fist exon ( the first array element )
  Vector_sort(models, Model_firstExonStartCompFunc);

  //my @models = sort { $a->[0]->start <=> $b->[0]->start }  @$models ;

  //# $model->[0]  = 1st exon 
  //# $model->[-3] = last exon 
  //# $model->[-2] = exon iuse
  //# $model->[-1] = score
  
  printf("Have %d to cluster\n", Vector_getNumElement(models));

  int startInd = 0;

  int i;
  for (i=0; i<Vector_getNumElement(models); i++) {
    Model *model = Vector_getElementAt(models, i);
    int clustered = 0;

    int nClust = Vector_getNumElement(clusters);
    Exon *firstExon = Vector_getElementAt(model->features, 0);
    Exon *lastExon  = Vector_getLastElement(model->features);

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

  printf("Have %d after clustering them\n", Vector_getNumElement(clusters));
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

  // the extra exon is a list of start end coords of the spliced intron sections
  // ie: end:start:end:start where the 1st and last coords are anchors to tie it 
  // into our rough model both must match before we can try and add any potentialy 
  // novel exons in
  Vector *sortedExons = Vector_copy(exons);
  Vector_sort(sortedExons, SeqFeature_startCompFunc);

  char **keys = StringHash_getKeys(extraExons);
  
  for (i=0; i<StringHash_getNumValues(extraExons); i++) {
    ExtraExonData *eed = StringHash_getValue(extraExons, keys[i]);

//    my @coords = split(/:/,$key);

//    long startAnchor = shift(@coords);
//    long endAnchor = pop(@coords);
// NOTE: I'm not shifting and popping (ie removing from array) so need to take account of this when making exons
    long startAnchor = eed->coords[0];
    long endAnchor   = eed->coords[eed->nCoord-1];

    //#print "START AND END $start_anchor  $end_anchor \n";

    if (startAnchor > endAnchor) {
      fprintf(stderr, "START AND END %ld %ld\n", startAnchor,  endAnchor);
      exit(1);
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
    if (RefineSolexaGenes_binSearchForOverlap(rsg, sortedExons, endAnchor)) {
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
  free(keys);

// NIY: Free sortedExons

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
        exit(1);
      }

      if (Exon_getEnd(rightExon) > Exon_getEnd(leftExon)) {
        Exon_setEnd(leftExon, Exon_getEnd(rightExon)); 
      }

      // get rid of right exon
      //splice(@exons,$i,1);
      Vector_removeElementAt(exons, i);

// NIY: Free removed exon??

      i--;
//#      @exons =  sort { $a->start <=> $b->start } @exons;
    }
  }
  
  printf("Merging exons - done overlap filter\n");
  
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
      if (strstr(DNAAlignFeature_getHitSeqName(intron), "non canonical")) continue;

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
      Vector_removeElementAt(exons, i-1);

// NIY: Free prevExon

      i--;
// Is this necessary???      next;
    }
  }
  printf("Merging exons done\n");
  
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
void RefineSolexaGenes_bamToIntronFeatures(RefineSolexaGenes *rsg, IntronBamConfig *intronBamConf, samfile_t *sam, bam_index_t *idx, int ref, int begRange, int endRange) {
  SliceAdaptor *sliceAdaptor = RefineSolexaGenes_getGeneSliceAdaptor(rsg);
  Vector *ifs = Vector_new();
  StringHash *extraExons = RefineSolexaGenes_getExtraExons(rsg);
  StringHash *idList = StringHash_new(STRINGHASH_SMALL);
  StringHash *readGroups = NULL;


  if (intronBamConf->groupNames != NULL && Vector_getNumElement(intronBamConf->groupNames) > 0) {
    
    printf("Limiting to read groups ");
    int i;
    readGroups = StringHash_new(STRINGHASH_SMALL);
    for (i=0; i<Vector_getNumElement(intronBamConf->groupNames); i++) {
      char *group = Vector_getElementAt(intronBamConf->groupNames, i);
      printf(" %s", group);
      StringHash_add(readGroups, group, &trueVal); 
    }
    printf("\n");
  }

  bam_iter_t iter = bam_iter_query(idx, ref, begRange, endRange);
  bam1_t *read = bam_init1();

  while (bam_iter_read(sam->x.bam, iter, read) >= 0) {

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

    Vector *mates = RefineSolexaGenes_getUngappedFeatures(rsg, read);

    Vector_sort(mates, CigarBlock_startCompFunc);
    //my @mates = sort { $a->[2] <=> $b->[2] } @{$self->ungapped_features($read)};

    // if mates > 2 then we have a possibility of adding in some extra exons into our rough models
    // as the read has spliced into and out of an exon
    // lets make them unique
    if (Vector_getNumElement(mates) > 2) {
      char keyString[2048]; keyString[0]='\0';
      long coords[1024]; // Hopefully won't have more than this!
      int nCoord = 0;

      int i;
      for (i=0; i<Vector_getNumElement(mates); i++) {
        CigarBlock *mate = Vector_getElementAt(mates, i);

        long start = mate->start;
        long end   = mate->end;

//  Unused      my $hstrand = $read->strand;

        if (i > 0) {
          sprintf(keyString, "%s%ld:", keyString, start);
          coords[nCoord++] = start;
        }
        if (i < Vector_getNumElement(mates)-1) {
          sprintf(keyString,"%s%ld:", keyString, end);
          coords[nCoord++] = end;
        }
      }

      if (!StringHash_contains(extraExons, keyString)) {
        StringHash_add(extraExons, keyString, ExtraExonData_new(coords, nCoord));
      }
      ExtraExonData *eed = StringHash_getValue(extraExons, keyString);

      eed->score++;
      //# print "Not doing extra_exon stuff for now\n";
    }

    //my $strand = $read->target->strand;
    int strand = bam1_strand(read) == 0 ? -1 : 1;
    if (intronBamConf->mixedBam) {
      if (spliced == '+') strand = 1;
      if (spliced == '-') strand = -1;
    } 

// Moved name setting out of loop
    char name[1024];
    strcpy(name, sam->header->target_name[read->core.tid]);
    StrUtil_strReplChr(name, '.', '*');

// Not used    long offset;
    char uniqueId[2048];
    int i;
    for (i=0; i<Vector_getNumElement(mates)-1; i++) {
      CigarBlock *mate   = Vector_getElementAt(mates, i);
      CigarBlock *mateP1 = Vector_getElementAt(mates, i+1);

      // intron reads should be split according to the CIGAR line
      // the default split function seems to ad
      // we want the ungapped features to make our introns
      // we dont allow . in the seq region name as we use them to delimit our paths
// Removed some unused vars here

      sprintf(uniqueId, "%s:%ld:%ld:%d", name, mate->end, mateP1->start, strand);

      if (!StringHash_contains(idList, uniqueId)) {
        StringHash_add(idList,
                       uniqueId,
                       IntronCoords_new(mate->end, mateP1->start, strand, -1, 0));
      }
      IntronCoords *ic = StringHash_getValue(idList, uniqueId);
      ic->score++;
    }
  }

/* SMJS These are for testing
  // For param testing store introns with different anal to results
  my $conslim = $ENV{CONSLIM};
  my $nonconslim = $ENV{NONCONSLIM};
  my $intron_anal = $self->create_analysis_object("intron_c" . $conslim . "_nc" . $nonconslim);
*/

  //# collapse them down and make them into simple features
  IntronCoords **icArray = StringHash_getValues(idList);
  Slice *chrSlice = RefineSolexaGenes_getChrSlice(rsg);
  Analysis *analysis = RefineSolexaGenes_getAnalysis(rsg);

  int i;
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
    if (length > 0) {
      sprintf(name,"%s:%ld:%ld:%d:", Slice_getName(chrSlice),
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

      if (!strcmp(Slice_getSeq(leftSplice), "NN") && !strcmp(Slice_getSeq(rightSplice), "NN")) {
        fprintf(stderr,"Warning: Cannot find dna sequence for %s this is used in detecting non cannonical splices\n", name);
      } else {
        //# is it cannonical
        if (DNAAlignFeature_getStrand(intFeat) == 1 ) {
          // is it GTAG?
          if (strcmp(Slice_getSeq(leftSplice), "GT") || strcmp(Slice_getSeq(rightSplice), "AG")) {
            canonical = 0;
          }
        } else {
          //# is it GTAG?
          if (strcmp(Slice_getSeq(rightSplice), "GT") || strcmp(Slice_getSeq(leftSplice), "AG")) {
            canonical = 0;
          }
        }
      }
      char hitName[2048];
      if (canonical) {
        sprintf(hitName, "%s%s", name, "canonical"); 
      } else {
        sprintf(hitName, "%s%s", name, "non canonical"); 
      }
      DNAAlignFeature_setHitSeqName(intFeat, hitName);
  
      Vector_addElement(ifs, intFeat);
// NIY: Free splice site sequence strings and slices??

    }
  }

  // sort them
  Vector_sort(ifs, SeqFeature_startCompFunc);
  if (RefineSolexaGenes_getFilterOnOverlapThreshold(rsg)) {
    Vector *tmpArray = Vector_new();
    int threshold = RefineSolexaGenes_getFilterOnOverlapThreshold(rsg);

    int arrayLength = Vector_getNumElement(ifs);
    if (arrayLength > 1) {
      int j;
      for (j=0; j<arrayLength-1; j++) {
        DNAAlignFeature *ifj = Vector_getElementAt(ifs, j);
        int k = 0;
        int count = 1;
        int overlappedSupport = 0;
        while (1) {
          ++k;
          DNAAlignFeature *ifjpk = Vector_getElementAt(ifs, j+k);

          if (count > threshold) {
            if (overlappedSupport < DNAAlignFeature_getScore(ifj)) {
//#            print STDERR "\t",$ifs[$j+$k]->hseqname, ': ', $ifs[$j+$k]->start, ':', $ifs[$j+$k]->end, "\n";
              Vector_addElement(tmpArray, ifj);
            } else {
//#          print STDERR 'THROWING: ', $ifs[$j]->hseqname, ': ', $ifs[$j]->start, ':', $ifs[$j]->end, "\n";
            }
            break;
          }
          overlappedSupport += DNAAlignFeature_getScore(ifjpk);

          if ((DNAAlignFeature_getEnd(ifj) < DNAAlignFeature_getStart(ifjpk)) || ((j+k) == arrayLength-1)) {
//#          print STDERR "\t",$ifs[$j+$k]->hseqname, ': ', $ifs[$j+$k]->start, ':', $ifs[$j+$k]->end, "\n";
            Vector_addElement(tmpArray, ifj);
            break;
          }
//#        print STDERR $ifs[$j+$k]->hseqname, "\n";
          if (DNAAlignFeature_getStrand(ifj) != DNAAlignFeature_getStrand(ifjpk)) {
            continue;
          }
          ++count;
        }
      }
// NIY: Free old ifs??
      ifs = tmpArray;
    }
  }
//#  print STDERR 'RES: ', scalar(@ifs), "\n";

/* Filtering code is Steve's experimental code
//# SMJS Filter here 

  if (!defined($conslim) || !defined($nonconslim)) {
    die "Env vars for CONSLIM and NONCONSLIM not set\n";
  }

  printf("Filter parameters:  Consensus splice coverage %f  Non consensus splice coverage %f\n", conslim, nonconslim);

  my @tmp;
  foreach my $f (@ifs) {
    if($f->hseqname =~ /non canonical/) {
      if ($f->score > $nonconslim && $f->end-$f->start < 50000) {
        push @tmp,$f;
      } else {
        print STDERR "Rejected non_canonical feature with score " . $f->score . "\n";
      }
    } else {
      if ($f->score > $conslim && $f->end-$f->start < 150000) {
        push @tmp,$f;
      } else {
        print STDERR "Rejected feature with score " . $f->score . "\n";
      }
    }
  }
  @ifs = @tmp;
*/

  RefineSolexaGenes_setIntronFeatures(rsg, ifs);

  RefineSolexaGenes_setExtraExons(rsg, extraExons);

  fprintf(stderr,"Got %d unique introns  ", Vector_getNumElement(ifs));
  fprintf(stderr," and %d potential novel exons from %s\n", StringHash_getNumValues(extraExons), intronBamConf->fileName);
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
Vector *RefineSolexaGenes_getUngappedFeatures(RefineSolexaGenes *rsg, bam1_t *b) {
  Vector *    ugfs = Vector_new();
  CigarBlock *lastMatchBlock;
  int         hadIntron = 0;
  int         nIntron = 0;
  CigarBlock *currentBlock = NULL;
  uint32_t *  cigar = bam1_cigar(b);

  int cigInd;
  int refPos;
  int readPos;
  for (cigInd = readPos = 0, refPos = b->core.pos+1; cigInd < b->core.n_cigar; ++cigInd) {
    int lenCigBlock = cigar[cigInd]>>4;
    int op          = cigar[cigInd]&0xf;

// M, =, X
    if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
      CigarBlock *block = CigarBlock_new(CB_MATCH, refPos, refPos+lenCigBlock-1);

      if (hadIntron) {
        hadIntron = 0;
        Vector_addElement(ugfs, block);
        currentBlock = NULL;
      } else {
        if (currentBlock) {
          CigarBlock_free(currentBlock);
        }
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
        Vector_addElement(ugfs, currentBlock);
        currentBlock = NULL;
      }
      refPos += lenCigBlock;
    } else {
      fprintf(stderr,"Cigar block type %d not supported\n", op);
      exit(1);
    }
  }
  if (hadIntron) {
    fprintf(stderr,"Error parsing cigar string - don't have two M regions surrounding an N region (intron)\n");
    exit(1);
  } 
  if (currentBlock) {
    CigarBlock_free(currentBlock);
  }
  if (Vector_getNumElement(ugfs) != nIntron+1) {
    fprintf(stderr,"Error parsing cigar string - don't have the expected number of blocks (%d) for %d introns (have %d)\n", 
            nIntron+1, nIntron, Vector_getNumElement(ugfs));
    exit(1);
  }

  return ugfs;

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
void RefineSolexaGenes_dnaToIntronFeatures(RefineSolexaGenes *rsg, long start, long end) {
// Unused  my $rough_genes = $self->prelim_genes;
  SliceAdaptor *intronSliceAdaptor = RefineSolexaGenes_getIntronSliceAdaptor(rsg);
  Slice *chrSlice = RefineSolexaGenes_getChrSlice(rsg);

  Slice *intronSlice = SliceAdaptor_fetchByRegion(intronSliceAdaptor, "toplevel", Slice_getSeqRegionName(chrSlice), start, end, 1, NULL, 0); 

  // featch all the dna_align_features for this slice by logic name
  Vector *logicNames = RefineSolexaGenes_getLogicNames(rsg);
  Vector *reads;

  if (Vector_getNumElement(logicNames)) {
    reads = Vector_new();
    fprintf(stderr,"Fetching reads with logic names: ");

    int i;
    for (i=0; i<Vector_getNumElement(logicNames); i++) {
      char *logicName = Vector_getElementAt(logicNames, i);

      fprintf(stderr,"%s ", logicName);
      Vector *feats = Slice_getAllDNAAlignFeatures(intronSlice, logicName, NULL, NULL, NULL);
      Vector_append(reads, feats); 
      Vector_free(feats);
    }
    fprintf(stderr,"\n");
  } else {
  //# fetch them all if no logic name is Supplied
    reads =  Slice_getAllDNAAlignFeatures(intronSlice, NULL, NULL, NULL, NULL);
  }
  fprintf(stderr, "Got %d reads\n", Vector_getNumElement(reads));

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
                DNAAlignFeature_getStart(ugfP1), DNAAlignFeature_getStrand(read), isCanonical ? "canonical" : "non canonical");

        if (!StringHash_contains(idList, uniqueId)) {
          StringHash_add(idList, 
                         uniqueId, 
                         IntronCoords_new(DNAAlignFeature_getEnd(ugf), DNAAlignFeature_getStart(ugfP1), DNAAlignFeature_getStrand(read), isCanonical, 0));
        }
        IntronCoords *ic = StringHash_getValue(idList, uniqueId);
        ic->score++;
      }
    }
  }

  fprintf(stderr, "Got %d collapsed introns\n", StringHash_getNumValues(idList));

  //# collapse them down and make them into dna align features


  Vector *intFeats = Vector_new();
  Analysis *analysis = RefineSolexaGenes_getAnalysis(rsg);

  IntronCoords **icArray = StringHash_getValues(idList);

  for (i=0; i<StringHash_getNumValues(idList); i++) {
    IntronCoords *ic = icArray[i];
    long length =  ic->nextExonStart - ic->prevExonEnd -1;

    char name[2048];
    if (length > 0) {
      sprintf(name,"%s:%ld:%ld:%d:%s", Slice_getName(chrSlice), 
                                       ic->prevExonEnd+1, 
                                       ic->nextExonStart-1, 
                                       ic->strand, 
                                       ic->isCanonical == 0 ? "non canonical" : "canonical");

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

      Vector_addElement(intFeats, intFeat);
    }
  }

  free(icArray);

  StringHash_free(idList, IntronCoords_free);

  // sort them
  Vector_sort(intFeats, SeqFeature_startCompFunc);

  RefineSolexaGenes_setIntronFeatures(rsg, intFeats);

  fprintf(stderr, "Got %d intron features\n", Vector_getNumElement(intFeats));
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

Vector *RefineSolexaGenes_fetchIntronFeatures(RefineSolexaGenes *rsg, long start, long end, long *offsetP) {
  Vector *sfs = RefineSolexaGenes_getIntronFeatures(rsg);

  long intronStart = 0;
  if (*offsetP) intronStart = *offsetP;
  int index = -1;

  Vector *chosenSf = Vector_new();

  // sfs is a sorted array
  int i;
  for (i=intronStart; i<Vector_getNumElement(sfs); i++) {
    DNAAlignFeature *intron = Vector_getElementAt(sfs, i);

    if (DNAAlignFeature_getStart(intron) > end) break;

    if (DNAAlignFeature_getStart(intron) <= end && DNAAlignFeature_getEnd(intron) >= start) {
      Vector_addElement(chosenSf, intron);

      // remember the position of the 1st intron to overlap
      // this exon - we will start counting from here next time
      if (index == -1) {
        index = i;
      }
    }
  }

  Vector *filteredIntrons = Vector_new();
// INTRON: 
  for (i=0; i<Vector_getNumElement(chosenSf); i++) {
    DNAAlignFeature *intron = Vector_getElementAt(chosenSf, i);

    if (strstr(DNAAlignFeature_getHitSeqName(intron), "non canonical")) {
      // check it has no overlap with any consensus introns
      // unless it out scores a consensus intron
      int done = 0;
      int j;
      for (j=0; j<Vector_getNumElement(chosenSf) && !done; j++) {
        DNAAlignFeature *compIntron = Vector_getElementAt(chosenSf, j);
        if (strstr(DNAAlignFeature_getHitSeqName(compIntron), "non canonical") == NULL) {
          if (DNAAlignFeature_getEnd(intron)    >  DNAAlignFeature_getStart(compIntron) && 
              DNAAlignFeature_getStart(intron)  <  DNAAlignFeature_getEnd(compIntron) &&
              DNAAlignFeature_getStrand(intron) == DNAAlignFeature_getStrand(compIntron)) {
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

  *offsetP = index;

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

  Exon *paddedExon = ExonUtils_createExon(start-20, end+20, -1, -1, -1, RefineSolexaGenes_getAnalysis(rsg), NULL, NULL, RefineSolexaGenes_getChrSlice(rsg));

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
  free(tmpVec);

  return paddedExon;
}

 
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
  if (rsg->intronFeatures != NULL) {
    fprintf(stderr, "Trying to set intronFeatures when its already set - may want append behaviour which isn't implemented - exiting\n");
    exit(1);
  }
// Perl did this sort
  Vector_sort(features, SeqFeature_startCompFunc);
  
  rsg->intronFeatures = features;
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

Vector *RefineSolexaGenes_getExtraExons(RefineSolexaGenes *rsg) {
  return rsg->extraExons;
}

void RefineSolexaGenes_setExtraExons(RefineSolexaGenes *rsg, Vector *extraExons) {
  rsg->extraExons = extraExons;
}


//####################################
//# config variable holders
//####################################

void RefineSolexaGenes_setIntronDb(RefineSolexaGenes *rsg, char *intronDb) {
  rsg->intronDb = StrUtil_copyString(&rsg->intronDb, intronDb, 0);
}

char *RefineSolexaGenes_getIntronDb(RefineSolexaGenes *rsg) {
  return rsg->intronDb;
}

void RefineSolexaGenes_setOutputDb(RefineSolexaGenes *rsg, char *outputDb) {
  rsg->outputDb = StrUtil_copyString(&rsg->outputDb, outputDb, 0);
}

char *RefineSolexaGenes_getOutputDb(RefineSolexaGenes *rsg) {
  return rsg->outputDb;
}

void RefineSolexaGenes_setModelDb(RefineSolexaGenes *rsg, char *modelDb) {
  rsg->modelDb = StrUtil_copyString(&rsg->modelDb, modelDb, 0);
}

char *RefineSolexaGenes_getModelDb(RefineSolexaGenes *rsg) {
  return rsg->modelDb;
}

void RefineSolexaGenes_setLogicName(RefineSolexaGenes *rsg, Vector *logicNames) {
  rsg->logicNames = logicNames;
}

Vector *RefineSolexaGenes_getLogicNames(RefineSolexaGenes *rsg) {
  return rsg->logicNames;
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
  rsg->bestScoreType = StrUtil_copyString(&rsg->bestScoreType, bestScoreType, 0);
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
  rsg->otherIsoformsType = StrUtil_copyString(&rsg->otherIsoformsType, otherIsoformsType, 0);
}

char *RefineSolexaGenes_getOtherIsoformsType(RefineSolexaGenes *rsg) {
  return rsg->otherIsoformsType;
}

void RefineSolexaGenes_setModelLogicName(RefineSolexaGenes *rsg, char *modelLN) {
  rsg->modelLogicName = StrUtil_copyString(&rsg->modelLogicName, modelLN, 0);
}

char *RefineSolexaGenes_getModelLogicName(RefineSolexaGenes *rsg) {
  return rsg->modelLogicName;
}

void RefineSolexaGenes_setBadModelsType(RefineSolexaGenes *rsg, char *badModelsType) {
  rsg->badModelsType = StrUtil_copyString(&rsg->badModelsType, badModelsType, 0);
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

void RefineSolexaGenes_setMinSingleExonCDSLength(RefineSolexaGenes *rsg, int minSingleExonCDSLength) {
  rsg->minSingleExonCDSLength = minSingleExonCDSLength;
}

int RefineSolexaGenes_getMinSingleExonCDSLength(RefineSolexaGenes *rsg) {
  return rsg->minSingleExonCDSLength;
}

void RefineSolexaGenes_setSingleExonModelType(RefineSolexaGenes *rsg, char *singleExonModelType) {
  rsg->singleExonModelType = StrUtil_copyString(&rsg->singleExonModelType, singleExonModelType, 0);
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

void RefineSolexaGenes_setWriteIntrons(RefineSolexaGenes *rsg, int writeIntrons) {
  rsg->writeIntrons = writeIntrons;
}

int RefineSolexaGenes_writeIntrons(RefineSolexaGenes *rsg) {
  return rsg->writeIntrons;
}

void RefineSolexaGenes_setTrimUtr(RefineSolexaGenes *rsg, int trimUtr) {
  rsg->trimUtr = trimUtr;
}

int RefineSolexaGenes_trimUtr(RefineSolexaGenes *rsg) {
  return rsg->trimUtr;
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

void RefineSolexaGenes_setRejectIntronCutoff(RefineSolexaGenes *rsg, double rejectIntronCutoff) {
  rsg->rejectIntronCutoff = rejectIntronCutoff;
}

double RefineSolexaGenes_getRejectIntronCutoff(RefineSolexaGenes *rsg) {
  return rsg->rejectIntronCutoff;
}

void RefineSolexaGenes_setDb(RefineSolexaGenes *rsg, DBAdaptor *db) {
  rsg->db = db;
}

DBAdaptor *RefineSolexaGenes_getDb(RefineSolexaGenes *rsg) {
  return rsg->db;
}


