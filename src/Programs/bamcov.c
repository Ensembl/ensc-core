/*
 * Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
 * Copyright [2016-2017] EMBL-European Bioinformatics Institute
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

/*
  Bamcov

  A program for calculating coverage values for features in a BAM format 

  Copyright (c) 1999-2013 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under an Apache 2 license.
  For license details, please see

    http://www.ensembl.org/info/about/legal/code_licence.html

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.
*/

#include <stdio.h>

#include "EnsC.h"

#include "DBAdaptor.h"
#include "BaseAdaptor.h"
#include "Basic/Vector.h"
#include "Slice.h"
#include "SliceAdaptor.h"
#include "StrUtil.h"
#include "IDHash.h"

#include "htslib/sam.h"
#include "htslib/hts.h"

void       Bamcov_usage();
int        calcCoverage(char *fName, Slice *slice, htsFile *in, hts_idx_t *idx, int flags);
int        geneStartCompFunc(const void *one, const void *two);
Vector *   getGenes(Slice *slice, int flags);

typedef struct CoverageStruct {
  long  coverage;
} Coverage;

// Flag values
#define M_UCSC_NAMING 1
#define M_VERIFY 4

// My custom bam core struct flag value
#define MY_FUSEDFLAG 32768

int verbosity = 1;

int main(int argc, char *argv[]) {
  DBAdaptor *      dba;
  StatementHandle *sth;
  ResultRow *      row;
  Vector *         slices;
  int              nSlices;
  htsFile *      out;

  int   argNum = 1;

  char *inFName  = NULL;
  char *outFName = NULL;

  char *dbUser = "ensro";
  char *dbPass = NULL;
  int   dbPort = 3306;

  char *dbHost = "ens-staging.internal.sanger.ac.uk";
  char *dbName = "homo_sapiens_core_71_37";

  char *assName = "GRCh37";

  char *chrName = "1";


  int flags = 0;
  int   threads  = 1;

  initEnsC(argc, argv);

  while (argNum < argc) {
    char *arg = argv[argNum];
    char *val;

// Ones without a val go here
    if (!strcmp(arg, "-U") || !strcmp(arg,"--ucsc_naming")) {
      flags |= M_UCSC_NAMING;
    } else {
// Ones with a val go in this block
      if (argNum == argc-1) {
        Bamcov_usage();
      }

      val = argv[++argNum];
  
      if (!strcmp(arg, "-i") || !strcmp(arg,"--in_file")) {
        StrUtil_copyString(&inFName,val,0);
      } else if (!strcmp(arg, "-o") || !strcmp(arg,"--out_file")) {
        StrUtil_copyString(&outFName,val,0);
      } else if (!strcmp(arg, "-h") || !strcmp(arg,"--host")) {
        StrUtil_copyString(&dbHost,val,0);
      } else if (!strcmp(arg, "-p") || !strcmp(arg,"--password")) {
        StrUtil_copyString(&dbPass,val,0);
      } else if (!strcmp(arg, "-P") || !strcmp(arg,"--port")) {
        dbPort = atoi(val);
      } else if (!strcmp(arg, "-n") || !strcmp(arg,"--name")) {
        StrUtil_copyString(&dbName,val,0);
      } else if (!strcmp(arg, "-u") || !strcmp(arg,"--user")) {
        StrUtil_copyString(&dbUser,val,0);
      } else if (!strcmp(arg, "-t") || !strcmp(arg,"--threads")) {
        threads = atoi(val);
      } else if (!strcmp(arg, "-a") || !strcmp(arg,"--assembly")) {
        StrUtil_copyString(&assName,val,0);
      } else if (!strcmp(arg, "-v") || !strcmp(arg,"--verbosity")) {
        verbosity = atoi(val);
// Temporary
      } else if (!strcmp(arg, "-c") || !strcmp(arg,"--chromosome")) {
        StrUtil_copyString(&chrName,val,0);
      } else {
        fprintf(stderr,"Error in command line at %s\n\n",arg);
        Bamcov_usage();
      }
    }
    argNum++;
  }

  if (verbosity > 0) {
    printf("Program for calculating read coverage in a BAM file \n"
           "Steve M.J. Searle.  searle@sanger.ac.uk  Last update April 2013.\n");
  }

  if (!inFName || !outFName) {
    Bamcov_usage();
  }

  dba = DBAdaptor_new(dbHost,dbUser,dbPass,dbName,dbPort,NULL);

  //nSlices = getSlices(dba, destName);
  nSlices = 1;

  slices = Vector_new();

  SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(dba);

  Slice *slice = SliceAdaptor_fetchByRegion(sa,NULL,chrName,POS_UNDEF,POS_UNDEF,1,NULL, 0);

  Vector_addElement(slices,slice);

  if (Vector_getNumElement(slices) == 0) {
    fprintf(stderr, "Error: No slices.\n");
    exit(1);
  }

  htsFile *in = hts_open(inFName, "rb");
  if (in == 0) {
    fprintf(stderr, "Fail to open BAM file %s\n", inFName);
    return 1;
  }

  hts_set_threads(in, threads);
  hts_idx_t *idx;
  idx = bam_index_load(inFName); // load BAM index
  if (idx == 0) {
    fprintf(stderr, "BAM index file is not available.\n");
    return 1;
  }

  int i;
  for (i=0; i<Vector_getNumElement(slices); i++) {
    Slice *slice = Vector_getElementAt(slices,i);

    if (verbosity > 0) printf("Working on '%s'\n",Slice_getName(slice));

//    if (verbosity > 0) printf("Stage 1 - retrieving annotation from database\n");
//    Vector *genes = getGenes(slice, flags);

    if (verbosity > 0) printf("Stage 1 - calculating coverage\n");
    calcCoverage(inFName, slice, in, idx, flags);
  }


  hts_idx_destroy(idx);
  hts_close(in);

  if (verbosity > 0) printf("Done\n");
  return 0;
}

/*
 Program usage message
*/
void Bamcov_usage() {
  printf("bamcov \n"
         "  -i --in_file     Input BAM file to map from (string)\n"
         "  -o --out_file    Output Wiggle file to write (string)\n"
         "  -U --ucsc_naming Input BAM file has 'chr' prefix on ALL seq region names (flag)\n"
         "  -h --host        Database host name for db containing genes (string)\n"
         "  -n --name        Database name for db containing genes (string)\n"
         "  -u --user        Database user (string)\n"
         "  -p --password    Database password (string)\n"
         "  -P --port        Database port (int)\n"
         "  -a --assembly    Assembly name (string)\n"
         "  -v --verbosity   Verbosity level (int)\n"
         "\n"
         "Notes:\n"
         "  -U will cause 'chr' to be prepended to all source seq_region names, except the special case MT which is changed to chrM.\n"
         "  -v Default verbosity level is 1. You can make it quieter by setting this to 0, or noisier by setting it > 1.\n"
         );
  exit(1);
}

int geneStartCompFunc(const void *one, const void *two) {
  Gene *g1 = *((Gene**)one);
  Gene *g2 = *((Gene**)two);

  return Gene_getStart(g1) - Gene_getStart(g2);
}

Vector *getGenes(Slice *slice, int flags) { 
  Vector *genes;
  genes = Slice_getAllGenes(slice, NULL, NULL, 1, NULL, NULL);

  Vector_sort(genes, geneStartCompFunc);

  return genes;
}

int calcCoverage(char *fName, Slice *slice, htsFile *in, hts_idx_t *idx, int flags) {
  int  ref;
  int  begRange;
  int  endRange;
  char region[1024];
  char region_name[512];


  if (Slice_getChrStart(slice) != 1) {
    fprintf(stderr, "Currently only allow a slice start position of 1\n");
    return 1;
  }
  if (flags & M_UCSC_NAMING) {
    sprintf(region,"chr%s", Slice_getSeqRegionName(slice));
  } else {
    sprintf(region,"%s", Slice_getSeqRegionName(slice));
  }
  bam_hdr_t *header = bam_hdr_init();
  header = bam_hdr_read(in->fp.bgzf);
  ref = bam_name2id(header, region);
  if (ref < 0) {
    fprintf(stderr, "Invalid region %s\n", region);
    exit(1);
  }
  sprintf(region,"%s:%ld-%ld", region_name,
                             Slice_getSeqRegionStart(slice),
                             Slice_getSeqRegionEnd(slice));
  if (hts_parse_reg(region, &begRange, &endRange) == NULL) {
    fprintf(stderr, "Could not parse %s\n", region);
    exit(2);
  }
  bam_hdr_destroy(header);


  hts_itr_t *iter = sam_itr_queryi(idx, ref, begRange, endRange);
  bam1_t *b = bam_init1();

  Coverage *coverage = calloc(Slice_getLength(slice),sizeof(Coverage));

  long counter = 0;
  long overlapping = 0;
  long bad = 0;
  int startIndex = 0;
  while (bam_itr_next(in, iter, b) >= 0) {
    if (b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP)) {
      bad++;
      continue;
    }

    int end;
    //end = bam_calend(&b->core, bam1_cigar(b));
    end = bam_endpos(b);

    // There is a special case for reads which have zero length and start at begRange (so end at begRange ie. before the first base we're interested in).
    // That is the reason for the || end == begRange test
    if (end == begRange) {
      continue;
    }
    counter++;

    if (!(counter%1000000)) {
      if (verbosity > 1) { printf("."); }
      fflush(stdout);
    }

// Remember: b->core.pos is zero based!
    int cigInd;
    int refPos;
    int readPos;
    uint32_t *cigar = bam_get_cigar(b);
    for (cigInd = readPos = 0, refPos = b->core.pos; cigInd < b->core.n_cigar; ++cigInd) {
      int k;
      int lenCigBlock = cigar[cigInd]>>4;
      int op          = cigar[cigInd]&0xf;

      if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
        for (k = 0; k < lenCigBlock; ++k) {
          //if (ref[refPos+k] == 0) break; // out of boundary
          coverage[refPos+k].coverage++;
        }
        if (k < lenCigBlock) break;
        refPos += lenCigBlock; readPos += lenCigBlock;
      } else if (op == BAM_CDEL) {
        for (k = 0; k < lenCigBlock; ++k) {
        //  if (ref[refPos+k] == 0) break;
          coverage[refPos+k].coverage++;
        }
        if (k < lenCigBlock) break;
        refPos += lenCigBlock;
      } else if (op == BAM_CSOFT_CLIP) {
        readPos += lenCigBlock;
      } else if (op == BAM_CHARD_CLIP) {
      } else if (op == BAM_CINS) {
         readPos += lenCigBlock;
      } else if (op == BAM_CREF_SKIP) {
         refPos += lenCigBlock;
      }
    }

#ifdef DONE
    int j;
    int done = 0;
    int hadOverlap = 0;
    
    for (j=startIndex; j < Vector_getNumElement(genes) && !done; j++) {
      Gene *gene = Vector_getElementAt(genes,j); 
      if (!gene) {
        continue;
      }
// Remember: b->core.pos is zero based!
      if (b->core.pos < Gene_getEnd(gene) && end >= Gene_getStart(gene)) {
        int k;

        int doneGene = 0;
        for (k=0; k<Gene_getTranscriptCount(gene) && !doneGene; k++) {
          Transcript *trans = Gene_getTranscriptAt(gene,k);

          if (b->core.pos < Transcript_getEnd(trans) && end >= Transcript_getStart(trans)) {
            int m;
     
            for (m=0; m<Transcript_getExonCount(trans) && !doneGene; m++) {
              Exon *exon = Transcript_getExonAt(trans,m);

              if (b->core.pos < Exon_getEnd(exon) && end >= Exon_getStart(exon)) {

                // Only count as overlapping once (could be that a read overlaps more than one gene)
                if (!hadOverlap) {
                  overlapping++;
                  hadOverlap = 1;
                }

                gs = IDHash_getValue(geneCountsHash, Gene_getDbID(gene));
                gs->score++;
                
                doneGene = 1;
              }
            }
          }
        }
      } else if (Gene_getStart(gene) > end) {
        done = 1;
      } else if (Gene_getEnd(gene) < b->core.pos+1) {
        gs = IDHash_getValue(geneCountsHash, Gene_getDbID(gene));
        printf("Gene %s (%s) score %ld\n",Gene_getStableId(gene), 
                                          Gene_getDisplayXref(gene) ? DBEntry_getDisplayId(Gene_getDisplayXref(gene)) : "", 
                                          gs->score);

        if (verbosity > 1) { 
          printf("Removing gene %s (index %d) with extent %d to %d\n", 
                 Gene_getStableId(gene), 
                 gs->index,
                 Gene_getStart(gene),
                 Gene_getEnd(gene));
        }
        Vector_setElementAt(genes,j,NULL);

        // Magic (very important for speed) - move startIndex to first non null gene
        int n;
        startIndex = 0;
        for (n=0;n<Vector_getNumElement(genes);n++) {
          void *v = Vector_getElementAt(genes,n);

          if (v != NULL) {
            break;
          }
          startIndex++;
        }
        if (verbosity > 1) { 
          printf("startIndex now %d\n",startIndex);
        }
      }
    }
#endif
  }
  if (verbosity > 1) { printf("\n"); }

#ifdef DONE
// Print out read counts for what ever's left in the genes array
  int n;
  for (n=0;n<Vector_getNumElement(genes);n++) {
    Gene *gene = Vector_getElementAt(genes,n);

    if (gene != NULL) {
      gs = IDHash_getValue(geneCountsHash, Gene_getDbID(gene));
      printf("Gene %s (%s) score %ld\n",Gene_getStableId(gene), 
                                        Gene_getDisplayXref(gene) ? DBEntry_getDisplayId(Gene_getDisplayXref(gene)) : "", 
                                        gs->score);
    }

  }
#endif

  printf("Read %ld reads. Number of bad reads (unmapped, qc fail, secondary, dup) %ld\n", counter, bad);

  long i;
  for (i=0; i< Slice_getLength(slice); i++) {
    printf("%ld %ld\n", i+1, coverage[i].coverage);
  }

  sam_itr_destroy(iter);
  bam_destroy1(b);


  return 1;
}

