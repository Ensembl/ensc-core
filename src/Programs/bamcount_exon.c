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
#include "DNAAlignFeature.h"
#include "StrUtil.h"
#include "IDHash.h"
#include "Transcript.h"

#include "sam.h"
#include "bam.h"

void       Bamcount_usage();
int        bamPosNameCompFunc(const void *one, const void *two);
int        bamPosCompFunc(const void *one, const void *two);
int        countReads(char *fName, Slice *slice, samfile_t *in, bam_index_t *idx, int flags, Vector *genes, IDHash *geneResultsHash, long long countUsableReads, FILE *bedFp);
long long  countReadsInFile(char *inFName);
bam1_t *   findMateInVector(bam1_t *b, Vector *vec);
int        findPosInVec(Vector *vec, int pos, char *bqname);
Vector *   flattenGene(Gene *gene);
int        geneStartCompFunc(const void *one, const void *two);
Vector *   getGenes(Slice *slice, int flags);
IDHash *   makeGeneResultsHash(Vector *genes);
bam1_t    *mateFoundInVectors(bam1_t *b, Vector **vectors);
void       printBam(FILE *fp, bam1_t *b, bam_header_t *header);

Vector *Bam_cigarToUngapped(bam1_t *b);

typedef struct GeneResultsStruct {
  int   index;
  long  score;
  Gene *gene;
  Vector *flatFeatures;
  Vector *exons;
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
  StatementHandle *sth;
  ResultRow *      row;
  Vector *         slices;
  int              nSlices;

  int   argNum = 1;

  char *inFName  = NULL;
  char *bedFName  = NULL;
  FILE *bedFp = NULL;

  char *dbUser = "ensro";
  char *dbPass = NULL;
  int   dbPort = 3306;

  char *dbHost = "ens-livemirror.internal.sanger.ac.uk";
  char *dbName = "homo_sapiens_core_71_37";

  char *assName = "GRCh37";

  char *chrName = "GL000206.1";


  int flags = 0;

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
      } else if (!strcmp(arg, "-b") || !strcmp(arg,"--bed")) {
        StrUtil_copyString(&bedFName,val,0);
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

  if (bedFName) {
    bedFp = fopen(bedFName, "w");
  }
  dba = DBAdaptor_new(dbHost,dbUser,dbPass,dbName,dbPort,NULL);

  //nSlices = getSlices(dba, destName);


  SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(dba);

  slices = SliceAdaptor_fetchAll(sa, "toplevel", NULL, 0);
  nSlices = Vector_getNumElement(slices);

//  Vector_setElementAt(slices, 0, Vector_getElementAt(slices, 26));
//  Vector_setNumElement(slices, 2);

//  slices = Vector_new();
//  Slice *slice = SliceAdaptor_fetchByRegion(sa,"chromosome", chrName,1,300000000,1, NULL, 0);
//  nSlices = 1;
//  Vector_addElement(slices,slice);
  

  if (Vector_getNumElement(slices) == 0) {
    fprintf(stderr, "Error: No slices.\n");
    exit(1);
  }


#ifdef _PBGZF_USE
  bam_set_num_threads_per(5);
#endif
//  long long totalUsableReads = countReadsInFile(inFName);
long long totalUsableReads = 1;

  fprintf(stderr, "Have %lld total usable reads\n",totalUsableReads);

#ifdef _PBGZF_USE
  bam_set_num_threads_per(5);
#endif
  samfile_t *in = samopen(inFName, "rb", 0);
  if (in == 0) {
    fprintf(stderr, "Fail to open BAM file %s\n", inFName);
    return 1;
  }

  bam_index_t *idx;
  idx = bam_index_load(inFName); // load BAM index
  if (idx == 0) {
    fprintf(stderr, "BAM index file is not available.\n");
    return 1;
  }

  int i;
  for (i=0; i<Vector_getNumElement(slices); i++) {
    Slice *slice = Vector_getElementAt(slices,i);

//    if (!strstr(Slice_getSeqRegionName(slice), "GL")) continue;
    if (verbosity > 0) fprintf(stderr, "Working on '%s'\n",Slice_getSeqRegionName(slice));

    if (verbosity > 0) fprintf(stderr, "Stage 1 - retrieving annotation from database\n");
    Vector *genes = getGenes(slice, flags);

    IDHash *geneResultsHash = makeGeneResultsHash(genes);

    if (verbosity > 0) fprintf(stderr, "Stage 2 - counting reads\n");
    countReads(inFName, slice, in, idx, flags, genes, geneResultsHash, totalUsableReads, bedFp);
  }


  bam_index_destroy(idx);
  samclose(in);

  if (verbosity > 0) fprintf(stderr, "Done\n");
  if (bedFp) fclose(bedFp);
  return 0;
}

long long countReadsInFile(char *inFName) {
  samfile_t *in = samopen(inFName, "rb", 0);
  long long nUsableReads = 0;

  if (in == 0) {
    fprintf(stderr, "Fail to open BAM file %s\n", inFName);
    exit(1);
  }

  bam1_t *b = bam_init1();

  int cnt = 0;
  while (bam_read1(in->x.bam, b) > 0 && b->core.tid >= 0) {
    if (!(b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP))) {
      nUsableReads++;
    }
    cnt++;
    if (!(cnt%1000000)) {
      fprintf(stderr,".");
      fflush(stdout);
    }
  }
  fprintf(stderr, "\n");

  samclose(in);

  return nUsableReads;
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
      return strcmp(bam1_qname(b1),bam1_qname(b2));
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

/*
 print out a bam1_t entry, particularly the flags (for debugging)
*/
void printBam(FILE *fp, bam1_t *b, bam_header_t *header) {    
  fprintf(fp, "%s %s %d %d %s (%d) %d %d %d\t\tP %d PP %d U %d MU %d R %d MR %d R1 %d R2 %d S %d QC %d D %d U %d\n",
                                  bam1_qname(b), 
                                  header->target_name[b->core.tid], 
                                  b->core.pos, 
                                  bam_calend(&b->core,bam1_cigar(b)),
                                  header->target_name[b->core.mtid], 
                                  b->core.mtid, 
                                  b->core.mpos, 
                                  b->core.isize, 
                                  bam_cigar2qlen(&(b->core),bam1_cigar(b)),
                                  b->core.flag & BAM_FPAIRED,
                                  b->core.flag & BAM_FPROPER_PAIR ? 1 : 0,
                                  b->core.flag & BAM_FUNMAP ? 1 : 0,
                                  b->core.flag & BAM_FMUNMAP ? 1 : 0,
                                  b->core.flag & BAM_FREVERSE ? 1 : 0,
                                  b->core.flag & BAM_FMREVERSE ? 1 : 0,
                                  b->core.flag & BAM_FREAD1 ? 1 : 0,
                                  b->core.flag & BAM_FREAD2 ? 1 : 0,
                                  b->core.flag & BAM_FSECONDARY ? 1 : 0,
                                  b->core.flag & BAM_FQCFAIL ? 1 : 0,
                                  b->core.flag & BAM_FDUP ? 1 : 0,
                                  b->core.flag & MY_FUSEDFLAG ? 1 : 0
                                  );
  fflush(fp);
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
    gr->exons = Gene_getAllExons(gene);
    
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
  Vector_free(exons);
 

  return blockFeatures;
}

int countReads(char *fName, Slice *slice, samfile_t *in, bam_index_t *idx, int flags, Vector *origGenesVec, IDHash *geneResultsHash, long long countUsableReads, FILE *bedFp) {
  int  ref;
  int  begRange;
  int  endRange;
  char region[1024];
  GeneResults *gr; 
  Vector *genes = Vector_copy(origGenesVec);
  //IDHash *transExonHash = IDHash_new(IDHASH_MEDIUM);

/*
  if (Slice_getSeqRegionStart(slice) != 1) {
    fprintf(stderr, "Currently only allow a slice start position of 1\n");
    return 1;
  }
*/
  if (flags & M_UCSC_NAMING) {
    sprintf(region,"chr%s:%ld-%ld", Slice_getSeqRegionName(slice), 
                                  Slice_getSeqRegionStart(slice), 
                                  Slice_getSeqRegionEnd(slice));
  } else {
    sprintf(region,"%s:%ld-%ld", Slice_getSeqRegionName(slice), 
                               Slice_getSeqRegionStart(slice), 
                               Slice_getSeqRegionEnd(slice));
  }
  bam_parse_region(in->header, region, &ref, &begRange, &endRange);
  if (ref < 0) {
    fprintf(stderr, "Invalid region %s %ld %ld\n", Slice_getSeqRegionName(slice), 
                                                 Slice_getSeqRegionStart(slice), 
                                                 Slice_getSeqRegionEnd(slice));
    return 1;
  }

  bam_iter_t iter = bam_iter_query(idx, ref, begRange, endRange);
  bam1_t *b = bam_init1();

  long counter = 0;
  long overlapping = 0;
  long bad = 0;
  int startIndex = 0;
  while (bam_iter_read(in->x.bam, iter, b) >= 0) {
    if (b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP)) {
      bad++;
      continue;
    }

    Vector *features = NULL;
    int end;
    end = bam_calend(&b->core, bam1_cigar(b));

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

        gr = IDHash_getValue(geneResultsHash, Gene_getDbID(gene));
        Vector *geneExons = gr->exons;
//        Vector *geneExons = Gene_getAllExons(gene);
        int m;
        for (m=0; m<Vector_getNumElement(geneExons) && !doneGene; m++) {
          Exon *exon = Vector_getElementAt(geneExons, m);

          if (b->core.pos < Exon_getEnd(exon) && end >= Exon_getStart(exon)) {
            if (features == NULL) {
              features = Bam_cigarToUngapped(b);
            }
            int n=0;
            for (n=0; n<Vector_getNumElement(features); n++) {
              DNAAlignFeature *daf = Vector_getElementAt(features, n);
              if (DNAAlignFeature_getStart(daf) <= Exon_getEnd(exon) && 
                  DNAAlignFeature_getEnd(daf) >= Exon_getStart(exon)) {
                Exon_setScore(exon, Exon_getScore(exon) + 1);
              }
            }
          }
        }
        //Vector_free(geneExons);

/*
        for (k=0; k<Gene_getTranscriptCount(gene) && !doneGene; k++) {
          Transcript *trans = Gene_getTranscriptAt(gene,k);
//          if (strcmp(Transcript_getBiotype(trans),"protein_coding")) {
//            continue;
//          }

          if (b->core.pos < Transcript_getEnd(trans) && end >= Transcript_getStart(trans)) {
            int m;
     
            for (m=0; m<Transcript_getExonCount(trans) && !doneGene; m++) {
              Exon *exon = Transcript_getExonAt(trans,m);

              if (!IDHash_contains(transExonHash, (IDType)exon)) {
                IDHash_add(transExonHash, (IDType) exon, trans);
              }
              if (IDHash_getValue(transExonHash, (IDType)exon) == trans) {
                if (b->core.pos < Exon_getEnd(exon) && end >= Exon_getStart(exon)) {
                  if (features == NULL) {
                    features = Bam_cigarToUngapped(b);
                  }
                  int n=0;
                  for (n=0; n<Vector_getNumElement(features); n++) {
                    DNAAlignFeature *daf = Vector_getElementAt(features, n);
                    if (DNAAlignFeature_getStart(daf) <= Exon_getEnd(exon) && 
                        DNAAlignFeature_getEnd(daf) >= Exon_getStart(exon)) {
                      Exon_setScore(exon, Exon_getScore(exon) + 1);
                    }
                  }
  
  //                gr = IDHash_getValue(geneResultsHash, Gene_getDbID(gene));
  //                gr->score++;
                  
                  //doneGene = 1;
                }
              }
            }
          }
        }
*/
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
    if (features) {
// Temporary - use free to free DNAAlignFeature
      Vector_setFreeFunc(features, free);
      Vector_free(features);
    }
  }
  if (verbosity > 1) { printf("\n"); }

// Print out read counts for what ever's left in the genes array
  int i;
  for (i=0;i<Vector_getNumElement(origGenesVec);i++) {
    Gene *gene = Vector_getElementAt(origGenesVec, i);

    gr = IDHash_getValue(geneResultsHash, Gene_getDbID(gene));

// Get max exon score
    double maxScore = 0.0;
    int j;
    for (j=0; j<Vector_getNumElement(gr->exons); j++) {
      Exon *exon = Vector_getElementAt(gr->exons, j);
      double normExonScore = Exon_getScore(exon)/Exon_getLength(exon);
      if (normExonScore > maxScore) {
        maxScore = normExonScore;
      }
    }
    fprintf(stderr, "Gene %s biotype %s num trans %d max normalised exon score %f\n", Gene_getStableId(gene),
            Gene_getBiotype(gene), Gene_getTranscriptCount(gene), maxScore);

    for (j=0; j<Gene_getTranscriptCount(gene); j++) {
      Transcript *trans = Gene_getTranscriptAt(gene, j);
      Transcript_sort(trans);
      fprintf(stderr, "Transcript %s biotype %s num exon %d\n", Transcript_getStableId(trans), 
              Transcript_getBiotype(trans), Transcript_getExonCount(trans));

      int k;
      double maxScoreInTrans = 0.0;
      for (k=0; k<Transcript_getExonCount(trans); k++) {
        Exon *exon = Transcript_getExonAt(trans, k);
        double normExonScore = Exon_getScore(exon)/Exon_getLength(exon);
        if (normExonScore > maxScoreInTrans) {
          maxScoreInTrans = normExonScore;
        }
      }

      Exon *prevExon = NULL;
      Exon *nextExon = NULL;
      for (k=0; k<Transcript_getExonCount(trans); k++) {
        Exon *exon = Transcript_getExonAt(trans, k);
        if (k<Transcript_getExonCount(trans)-1) {
          nextExon = Transcript_getExonAt(trans, k+1);
        } else {
          nextExon = NULL;
        }

        double normExonScore = Exon_getScore(exon)/Exon_getLength(exon);
        fprintf(stderr, "Exon score for %s (%s) %s gsource %s %s %s tsource %s %s %s %ld %ld %d length %ld score %f score per base %f norm score/base %f", 
                Gene_getStableId(gene), 
                Gene_getDisplayXref(gene) ? DBEntry_getDisplayId(Gene_getDisplayXref(gene)) : "", 
                Gene_getBiotype(gene), Gene_getSource(gene),
                Transcript_getStableId(trans), Transcript_getBiotype(trans), Analysis_getLogicName(Transcript_getAnalysis(trans)), 
                Exon_getStableId(exon), Slice_getSeqRegionName(Exon_getSlice(exon)), Exon_getSeqRegionStart(exon), Exon_getSeqRegionEnd(exon), Exon_getSeqRegionStrand(exon), 
                Exon_getLength(exon), Exon_getScore(exon), Exon_getScore(exon)/Exon_getLength(exon), 
                maxScore != 0.0 ? (((normExonScore)/maxScore)*1000.0) : 0.0);

        if (normExonScore < maxScoreInTrans/15.0 && prevExon!=NULL && nextExon!=NULL && Exon_getLength(exon) > 50) {
          double prevNormExonScore = Exon_getScore(prevExon)/Exon_getLength(prevExon);
          double nextNormExonScore = Exon_getScore(nextExon)/Exon_getLength(nextExon);
          if (prevNormExonScore/5.0 > normExonScore && nextNormExonScore/5.0 > normExonScore) {
            fprintf(stderr," LOWSCORE");
          }
        }
        fprintf(stderr,"\n");

        prevExon = exon;
      }
    }

/*
// Not really rpkm because only reads from this chr
    gr = IDHash_getValue(geneResultsHash, Gene_getDbID(gene));

    double rpkm = gr->score * 1000000000.0 / countUsableReads;
*/

/*
    printf("Gene %s (%s) score %ld flatlength %ld rpkm (sort of) %-14.5f\n",Gene_getStableId(gene), 
                                        Gene_getDisplayXref(gene) ? DBEntry_getDisplayId(Gene_getDisplayXref(gene)) : "", 
                                        gr->score, gr->flatLength, rpkm);
*/
/*
    printf("Gene %s (%s) score %ld flatlength %ld\n",Gene_getStableId(gene), 
                                        Gene_getDisplayXref(gene) ? DBEntry_getDisplayId(Gene_getDisplayXref(gene)) : "", 
                                        gr->score, gr->flatLength);
*/
  }

/*
  printf("Read %ld reads. Num overlapping exons %ld. Number of bad reads (unmapped, qc fail, secondary, dup) %ld\n", counter, overlapping, bad);
*/

  bam_iter_destroy(iter);
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
  char *bqname = bam1_qname(b);
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
            
          int qnameCmp =  strcmp(bam1_qname(vb),bqname); // Note order of comparison
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
      //int qnameCmp =  strcmp(bqname,bam1_qname(b)); // Note order of comparison
      //if (!qnameCmp) {
        // key found at index imid
        // Back up through array to find first index which matches (can be several)
        for (;imid>=0;imid--) {
          bam1_t *vb = Vector_getElementAt(vec,imid);
       //   if (vb->core.pos != pos || strcmp(bqname,bam1_qname(vb))) {
          if (vb->core.pos != pos) {// || strcmp(bqname,bam1_qname(vb))) {
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

Vector *Bam_cigarToUngapped(bam1_t *b) {
  int cigInd;
  int refPos;
  int readPos;
  uint32_t *cigar = bam1_cigar(b);

 Vector *features = Vector_new();

  for (cigInd = readPos = 0, refPos = b->core.pos; cigInd < b->core.n_cigar; ++cigInd) {
    int k;
    int lenCigBlock = cigar[cigInd]>>4;
    int op          = cigar[cigInd]&0xf;

    if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
      DNAAlignFeature *daf = DNAAlignFeature_new();

      DNAAlignFeature_setStart(daf, refPos+1);
      DNAAlignFeature_setEnd(daf, refPos+lenCigBlock+1);

      Vector_addElement(features, daf);

      refPos += lenCigBlock; readPos += lenCigBlock;
    } else if (op == BAM_CDEL) {
/* Not sure if I want a feature for this
      for (k = 0; k < lenCigBlock; ++k) {
      //  if (ref[refPos+k] == 0) break;
        coverage[refPos+k].coverage++;
      }
      if (k < lenCigBlock) break;
*/
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

  return features;
}
