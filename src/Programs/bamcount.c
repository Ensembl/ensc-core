/*
 * Copyright [1999-2024] EMBL-European Bioinformatics Institute
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
  Bamcount

  A program for calculating RPKM type values for features in a BAM format 
  against a set of annotation from an Ensembl database

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
#include "SliceAdaptor.h"
#include "Slice.h"
#include "StrUtil.h"
#include "IDHash.h"
#include "Transcript.h"

#include "bamhelper.h"
#include "htslib/sam.h"
#include "htslib/hts.h"

void       Bamcount_usage();
int        bamPosNameCompFunc(const void *one, const void *two);
int        bamPosCompFunc(const void *one, const void *two);
int        countReads(char *fName, Slice *slice, htsFile *in, hts_idx_t *idx, int flags, Vector *genes, IDHash *geneResultsHash, long long countUsableReads);
bam1_t *   findMateInVector(bam1_t *b, Vector *vec);
int        findPosInVec(Vector *vec, int pos, char *bqname);
Vector *   flattenGene(Gene *gene);
int        geneStartCompFunc(const void *one, const void *two);
Vector *   getGenes(Slice *slice, int flags);
IDHash *   makeGeneResultsHash(Vector *genes);
bam1_t    *mateFoundInVectors(bam1_t *b, Vector **vectors);

typedef struct GeneResultsStruct {
  int   index;
  long  score;
  Gene *gene;
  Vector *flatFeatures;
  long  flatLength;
} GeneResults;

// Flag values
#define M_UCSC_NAMING 1
#define M_VERIFY 4

// My custom bam core struct flag value
#define MY_FUSEDFLAG 32768

int verbosity = 1;

int main(int argc, char *argv[]) {
  DBAdaptor *      dba;
  Vector *         slices;
  int              nSlices;

  int   argNum = 1;

  char *inFName  = NULL;

  char *dbUser = "ensro";
  char *dbPass = NULL;
  int   dbPort = 3306;

  char *dbHost = "ens-livemirror.internal.sanger.ac.uk";
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
        Bamcount_usage();
      }

      val = argv[++argNum];
  
      if (!strcmp(arg, "-i") || !strcmp(arg,"--in_file")) {
        StrUtil_copyString(&inFName,val,0);
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
        Bamcount_usage();
      }
    }
    argNum++;
  }

  if (verbosity > 0) {
    printf("Program for read counting of BAM reads against annotation \n"
           "Steve M.J. Searle.  searle@sanger.ac.uk  Last update Mar 2013.\n");
  }

  if (!inFName) {
    Bamcount_usage();
  }

  dba = DBAdaptor_new(dbHost,dbUser,dbPass,dbName,dbPort,NULL);

  //nSlices = getSlices(dba, destName);
  nSlices = 1;

  slices = Vector_new();

  SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(dba);

  Slice *slice = SliceAdaptor_fetchByRegion(sa,"chromosome", chrName,1,300000000,1, NULL, 0);

  Vector_addElement(slices,slice);

  if (Vector_getNumElement(slices) == 0) {
    fprintf(stderr, "Error: No slices.\n");
    exit(1);
  }


//  long long totalUsableReads = countReadsInFile(inFName);
long long totalUsableReads = 1;

  printf("Have %lld total usable reads\n",totalUsableReads);

  htsFile *in = hts_open(inFName, "rb");
  if (in == 0) {
    fprintf(stderr, "Fail to open BAM file %s\n", inFName);
    return 1;
  }

  hts_set_threads(in, threads);
  hts_idx_t *idx;
  idx = sam_index_load(in, inFName); // load BAM index
  if (idx == 0) {
    fprintf(stderr, "BAM index file is not available.\n");
    return 1;
  }

  int i;
  for (i=0; i<Vector_getNumElement(slices); i++) {
    Slice *slice = Vector_getElementAt(slices,i);

    if (verbosity > 0) printf("Working on '%s'\n",Slice_getSeqRegionName(slice));

    if (verbosity > 0) printf("Stage 1 - retrieving annotation from database\n");
    Vector *genes = getGenes(slice, flags);

    IDHash *geneResultsHash = makeGeneResultsHash(genes);

    if (verbosity > 0) printf("Stage 2 - counting reads\n");
    countReads(inFName, slice, in, idx, flags, genes, geneResultsHash, totalUsableReads);
  }


  hts_idx_destroy(idx);
  hts_close(in);

  if (verbosity > 0) printf("Done\n");
  return 0;
}

/*
 Program usage message
*/
void Bamcount_usage() {
  printf("bamcount \n"
         "  -i --in_file     Input BAM file to map from (string)\n"
         "  -U --ucsc_naming Input BAM file has 'chr' prefix on ALL seq region names (flag)\n"
         "  -h --host        Database host name for db containing mapping (string)\n"
         "  -n --name        Database name for db containing mapping (string)\n"
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

int bamPosNameCompFunc(const void *one, const void *two) {
  bam1_t *b1 = *((bam1_t**)one);
  bam1_t *b2 = *((bam1_t**)two);

  if (b1->core.tid < b2->core.tid) {
    printf("Urrr.... this shouldn't happen b1 tid %d   b2 tid %d\n", b1->core.tid, b2->core.tid);
    return -1;
  } else if (b1->core.tid > b2->core.tid) {
    printf("Urrr.... this shouldn't happen b1 tid %d   b2 tid %d\n", b1->core.tid, b2->core.tid);
    return 1;
  } else {
    if (b1->core.pos > b2->core.pos) {
      return 1;
    } else if (b1->core.pos < b2->core.pos) {
      return -1;
    } else {
      return strcmp(bam_get_qname(b1),bam_get_qname(b2));
    }
  }
}


int bamPosCompFunc(const void *one, const void *two) {
  bam1_t *b1 = *((bam1_t**)one);
  bam1_t *b2 = *((bam1_t**)two);

  return b1->core.pos - b2->core.pos;
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

// Pregenerate the hash for gene scores, so we don't have to test for existence within the BAM read loop
// Also means all genes will have a score hash entry
IDHash *makeGeneResultsHash(Vector *genes) {
  IDHash *geneResultsHash = IDHash_new(IDHASH_LARGE);

  int i;
  for (i=0; i < Vector_getNumElement(genes); i++) {
    Gene *gene = Vector_getElementAt(genes,i);
    GeneResults *gr;

    if ((gr = (GeneResults *)calloc(1,sizeof(GeneResults))) == NULL) {
      fprintf(stderr,"ERROR: Failed allocating GeneResults\n");
      exit(1);
    }
    gr->index = i;
    gr->gene  = gene;
    gr->flatFeatures = flattenGene(gene);
    int j=0;
    gr->flatLength = 0;
    for (j=0; j<Vector_getNumElement(gr->flatFeatures); j++) {
      SeqFeature *sf = Vector_getElementAt(gr->flatFeatures,j);
      gr->flatLength += SeqFeature_getLength(sf); 
    }
    
    IDHash_add(geneResultsHash,Gene_getDbID(gene),gr);
  }

  return geneResultsHash;
}

Vector *flattenGene(Gene *gene) {
  SeqFeature *curBlock = NULL;
  Vector *exons = Gene_getAllExons(gene);
  Vector *blockFeatures = Vector_new();

  Vector_sort(exons, SeqFeature_startCompFunc);

  int i;
  for (i=0; i<Vector_getNumElement(exons); i++) {
    Exon *exon = Vector_getElementAt(exons, i);
    long exStart = Exon_getStart(exon) <= 0 ? 1 : Exon_getStart(exon);
    long exEnd   = Exon_getEnd(exon) <= 0 ?   1 : Exon_getEnd(exon);

    if (curBlock && SeqFeature_getEnd(curBlock) >= exStart) {
      // printf("Adding to block with range %ld to %ld\n", exStart, exEnd);
      // printf("          Block was %d to %d\n", SeqFeature_getStart(curBlock), SeqFeature_getEnd(curBlock));
      if (exEnd > SeqFeature_getEnd(curBlock)) SeqFeature_setEnd(curBlock, exEnd);
    } else {
      // printf("Starting new block with range %ld to %ld\n", exStart, exEnd);
      curBlock = SeqFeature_new();
      SeqFeature_setStart(curBlock, exStart);
      SeqFeature_setEnd(curBlock, exEnd);
      Vector_addElement(blockFeatures, curBlock);
    }
  }

  Vector_sort(blockFeatures, SeqFeature_startCompFunc);

  return blockFeatures;
}

int countReads(char *fName, Slice *slice, htsFile *in, hts_idx_t *idx, int flags, Vector *origGenesVec, IDHash *geneResultsHash, long long countUsableReads) {
  int  ref;
  int  begRange;
  int  endRange;
  char region[1024];
  char region_name[512];
  GeneResults *gr; 
  Vector *genes = Vector_copy(origGenesVec);


  if (Slice_getSeqRegionStart(slice) != 1) {
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
          if (strcmp(Transcript_getBiotype(trans),"protein_coding")) {
            continue;
          }

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

                gr = IDHash_getValue(geneResultsHash, Gene_getDbID(gene));
                gr->score++;
                
                doneGene = 1;
              }
            }
          }
        }
      } else if (Gene_getStart(gene) > end) {
        done = 1;
      } else if (Gene_getEnd(gene) < b->core.pos+1) {
//        gr = IDHash_getValue(geneResultsHash, Gene_getDbID(gene));
//        printf("Gene %s (%s) score %ld\n",Gene_getStableId(gene), 
//                                          Gene_getDisplayXref(gene) ? DBEntry_getDisplayId(Gene_getDisplayXref(gene)) : "", 
//                                          gr->score);

        if (verbosity > 1) { 
          printf("Removing gene %s (index %d) with extent %ld to %ld\n", 
                 Gene_getStableId(gene), 
                 gr->index,
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
  }
  if (verbosity > 1) { printf("\n"); }

// Print out read counts for what ever's left in the genes array
  int n;
  for (n=0;n<Vector_getNumElement(origGenesVec);n++) {
    Gene *gene = Vector_getElementAt(origGenesVec, n);
// Not really rpkm because only reads from this chr

    gr = IDHash_getValue(geneResultsHash, Gene_getDbID(gene));

    double rpkm = gr->score * 1000000000.0 / countUsableReads;

    printf("Gene %s (%s) score %ld flatlength %ld rpkm (sort of) %-14.5f\n",Gene_getStableId(gene), 
                                        Gene_getDisplayXref(gene) ? DBEntry_getDisplayId(Gene_getDisplayXref(gene)) : "", 
                                        gr->score, gr->flatLength, rpkm);
  }

  printf("Read %ld reads. Num overlapping exons %ld. Number of bad reads (unmapped, qc fail, secondary, dup) %ld\n", counter, overlapping, bad);

  sam_itr_destroy(iter);
  bam_destroy1(b);
}

bam1_t *getMateFromRemoteMates(bam1_t *b, Vector **remoteMates) {
// For now just delegate, could optimise more
  return mateFoundInVectors(b, remoteMates);
}

bam1_t *mateFoundInVectors(bam1_t *b, Vector **vectors) {
  Vector *vec = vectors[b->core.mtid];
  
  if (!vec) {
    return NULL;
  }

  return findMateInVector(b, vec);
}

bam1_t *findMateInVector(bam1_t *b, Vector *vec) {
  char *bqname = bam_get_qname(b);
  int firstInd = findPosInVec(vec, b->core.mpos, bqname);
  int i;

  if (firstInd >=0 ) {
    for (i=firstInd;i<Vector_getNumElement(vec);i++) {
      bam1_t *vb = (bam1_t *)Vector_getElementAt(vec,i);

      if (b->core.mpos == vb->core.pos) {
        //if (b->core.isize == -vb->core.isize &&
        if (abs(b->core.isize) == abs(vb->core.isize) &&
            b->core.pos == vb->core.mpos &&
            !(vb->core.flag & MY_FUSEDFLAG)) {
            
          int qnameCmp =  strcmp(bam_get_qname(vb),bqname); // Note order of comparison
          if (!qnameCmp) { // Name match
            return vb;
          } else if (qnameCmp > 0) { // Optimisation: As names in vector are sorted, once string comparison is positive can not be any matches anymore
            break;
          }
        }
      } else if (b->core.mpos < vb->core.pos) {
        break;
      }
    }
  }

  return NULL;
}

int findPosInVec(Vector *vec, int pos, char *bqname) {
  int imin = 0;
  int imax = Vector_getNumElement(vec)-1;

  while (imax >= imin) {
    int imid = (imax+imin) / 2;
    bam1_t *b = Vector_getElementAt(vec,imid);

//    printf("imid = %d imin = %d imax = %d pos = %d b->core.pos = %d\n",imid,imin,imax,pos,b->core.pos);

    if (pos > b->core.pos) {
      imin = imid + 1;
    } else if (pos < b->core.pos) {
      imax = imid - 1;
    } else {
      //int qnameCmp =  strcmp(bqname,bam_get_qname(b)); // Note order of comparison
      //if (!qnameCmp) {
        // key found at index imid
        // Back up through array to find first index which matches (can be several)
        for (;imid>=0;imid--) {
          bam1_t *vb = Vector_getElementAt(vec,imid);
       //   if (vb->core.pos != pos || strcmp(bqname,bam_get_qname(vb))) {
          if (vb->core.pos != pos) {// || strcmp(bqname,bam_get_qname(vb))) {
            return imid+1;
          }
        }
      //} else if (qnameCmp < 0) {
      //  imax = imid - 1;
      //} else {
      //  imin = imid + 1;
      //}

      return 0; // Is the first element
    }
  }
  return -1;
}

