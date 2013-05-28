/*
  Bammap

  A program for mapping features in a BAM format file between to assemblies.
  Input is a BAM file with features on the source coordinates
  Output is a BAM file with features on the destination coordinates
  It uses mapping data stored in an Ensembl database.

  It is a fairly efficient program (reads, maps and writes around 250 million 
  reads per hour on my test system). 

  Although I've tried to limit memory use as much as possible, it can require 
  fairly large amounts of memory, particularly if there are many pairs on the 
  source coordinate system which map to different mapping blocks (remote mates),
  or where many reads only partially map to a mapping block. The largest file I've
  mapped is around 200Gb, which used around 21.3Gb of RAM whilst running.

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
#include "CoordSystemAdaptor.h"
#include "BaseAdaptor.h"
#include "Basic/Vector.h"
#include "Slice.h"
#include "StrUtil.h"

#include "sam.h"
#include "bam.h"

typedef struct mappingStruct {
  Slice *sourceSlice;
  Slice *destSlice;
  int    ori;
} Mapping;

typedef struct readMapStatsStruct {
  long long nRead;
  long long nWritten;
  long long nOverEnds;
  long long nRemoteMate;
  long long nReversed;
  long long nUnmappedMate;
} ReadMapStats;

void       Bammap_usage();
int        bamPosNameCompFunc(const void *one, const void *two);
int        bamPosCompFunc(const void *one, const void *two);
int        calcNewEnd(Mapping *mapping, bam1_t *b, int sourceEnd);
int        calcNewPos(Mapping *mapping, bam1_t *b, int sourceEnd);
void       clearPairing(bam1_t *b, int end);
bam1_t *   findMateInVector(bam1_t *b, Vector *vec);
int        findPosInVec(Vector *vec, int pos, char *bqname);
Vector *   getDestinationSlices(DBAdaptor *dba, char *assName);
Vector *   getMappings(DBAdaptor *dba, char *seqName, char *fromAssName, char *toAssName, int rev, int flags);
Vector **  getMappingVectorsBySourceRegion(DBAdaptor *dba, samfile_t *in, char *sourceName, char *destName, int flags);
bam1_t *   getMateFromRemoteMates(bam1_t *b, Vector **remoteMates);
int        getPairedMappingFailLists(samfile_t *in, Vector **mappingVectors, Vector ***failedVectorsP, Vector ***remoteMatesP);
int        mapBam(char *fName, samfile_t *out, Mapping *mapping, ReadMapStats *regionStats, samfile_t *in, 
                  bam_index_t *idx, Vector **mappingVectors, Vector **failedVectors, Vector **remoteMates, int flags);
int        mapBam_forward(Mapping *mapping, bam1_t *b, samfile_t *in, samfile_t *out, Vector **mappingVectors, Vector **failedVectors, Vector **remoteMates, 
                          int begRange, int endRange, ReadMapStats *regionStats, int newtid, bam_iter_t iter);
int        mapBam_reverse(Mapping *mapping, bam1_t *b, samfile_t *in, samfile_t *out, Vector **mappingVectors, Vector **failedVectors, Vector **remoteMates, 
                          int begRange, int endRange, ReadMapStats *regionStats, int newtid, bam_iter_t iter);
int        mapLocation(Mapping *mapping, int pos);
int        mapMateLocation(Mapping *mapping, bam1_t *b, int end, samfile_t *in, samfile_t *out, Vector **mappingVectors, 
                           Vector **failedVectors, Vector **remoteMates, Vector *reverseCache, 
                           ReadMapStats *regionStats, int begRange, int endRange, int newtid, int newpos);
int        mapRemoteLocation(Vector **mappingVectors, int seqid, int pos, Mapping **containingMappingP);
bam1_t    *mateFoundInVectors(bam1_t *b, Vector **vectors);
void       printBam(FILE *fp, bam1_t *b, bam_header_t *header);
void       printMapping(FILE *fp, Mapping *mapping);
int        verifyBam(samfile_t *sam);
samfile_t *writeBamHeader(char *inFName, char *outFName, Vector *destinationSlices);

// Flag values
#define M_UCSC_NAMING 1
//#define M_QUICK 2
#define M_VERIFY 4

// My custom bam core struct flag value
#define MY_FUSEDFLAG 32768

int verbosity = 1;

int main(int argc, char *argv[]) {
  DBAdaptor *      dba;
  StatementHandle *sth;
  ResultRow *      row;
  Vector *         destinationSlices;
  Vector *         mappings;
  samfile_t *      out;

  int   argNum = 1;

  char *inFName  = NULL;
  char *outFName = NULL;

  char *dbUser = "ensro";
  char *dbPass = NULL;
  int   dbPort = 3306;

  char *dbHost = "ens-livemirror.internal.sanger.ac.uk";
  char *dbName = "homo_sapiens_core_70_37";

  char *sourceName = "GRCh37";
  char *destName   = "NCBI36";

//  char *dbHost = "genebuild2.internal.sanger.ac.uk";
//  char *dbName = "ba1_nomascus_leucogenys_core_70_1_nleu3";

//  char *sourceName = "Nleu1.0";
//  char *destName   = "Nleu_3.0";

  int flags = 0;

  ReadMapStats totalStats;
  ReadMapStats regionStats;


  memset(&totalStats, 0, sizeof(ReadMapStats));

  initEnsC(argc, argv);

  while (argNum < argc) {
    char *arg = argv[argNum];
    char *val;

// Ones without a val go here
// Disable -q because complicates mapping logic and doesn't speed it up much. If need
// to implement it (maybe to save memory), mapBam remote and local reverse mappings need changing
//    if (!strcmp(arg, "-q") || !strcmp(arg,"--quick")) {
//      flags |= M_QUICK;
//    } else 
    if (!strcmp(arg, "-U") || !strcmp(arg,"--ucsc_naming")) {
      flags |= M_UCSC_NAMING;
    } else if (!strcmp(arg, "-V") || !strcmp(arg,"--verify")) {
      flags |= M_VERIFY;
    } else {
// Ones with a val go in this block
      if (argNum == argc-1) {
        Bammap_usage();
      }
  
      val = argv[++argNum];
  //    printf("%s %s\n",arg,val);
  
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
      } else if (!strcmp(arg, "-s") || !strcmp(arg,"--source_ass")) {
        StrUtil_copyString(&sourceName,val,0);
      } else if (!strcmp(arg, "-d") || !strcmp(arg,"--dest_ass")) {
        StrUtil_copyString(&destName,val,0);
      } else if (!strcmp(arg, "-v") || !strcmp(arg,"--verbosity")) {
        verbosity = atoi(val);
      } else {
        fprintf(stderr,"Error in command line at %s\n\n",arg);
        Bammap_usage();
      }
    }

    argNum++;
  }

  if (verbosity > 0) {
    printf("Program for mapping BAM files between assemblies\n"
           "Steve M.J. Searle.  searle@sanger.ac.uk  Last update Jan 2013.\n");
  }

  if (!inFName || !outFName) {
    Bammap_usage();
  }

  if (verbosity > 0) printf("Stage 1 - retrieving %s to %s mappings from database\n",sourceName,destName);
  dba = DBAdaptor_new(dbHost,dbUser,dbPass,dbName,dbPort,NULL);

  destinationSlices = getDestinationSlices(dba, destName);

  if (Vector_getNumElement(destinationSlices) == 0) {
    fprintf(stderr, "Error: No destination slices.\n");
    exit(1);
  }

  out = writeBamHeader(inFName,outFName,destinationSlices);

//  bam_flush(out->x.bam);

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

  // Load all mappings into an array of vectors arranged by source region id
  // Only used for remote mate mapping, where we don't know in advanced which
  // region the mapping will be to
  Vector **mappingVectors = getMappingVectorsBySourceRegion(dba, in, sourceName, destName, flags);

  if (verbosity > 0) printf("Stage 2 - search for reads which don't map\n");

  Vector **failedVectors;
  Vector **remoteMates;
  //if (!(flags & M_QUICK)) {
  getPairedMappingFailLists(in, mappingVectors, &failedVectors, &remoteMates);
  //} else {
  //  failedVectors = calloc(in->header->n_targets, sizeof(Vector *));
  //  remoteMates   = calloc(in->header->n_targets, sizeof(Vector *));
  //}


  if (verbosity > 0) printf("Stage 3 - reading, transforming and writing all mapping reads\n");
  int i;
  for (i=0; i<Vector_getNumElement(destinationSlices); i++) {
    Slice *slice = Vector_getElementAt(destinationSlices,i);

    if (verbosity > 0) printf("Working on '%s'\n",Slice_getChrName(slice));
    //if (!strcmp(Slice_getChrName(slice),"3")) break;

    mappings = getMappings(dba,Slice_getChrName(slice),sourceName,destName,0, flags);
    int j;
    for (j=0;j<Vector_getNumElement(mappings); j++) {
      Mapping *mapping = Vector_getElementAt(mappings,j);
      if (verbosity > 1) printMapping(stdout, mapping);
      
      memset(&regionStats, 0, sizeof(ReadMapStats));

      mapBam(inFName, out, mapping, &regionStats, in, idx, mappingVectors, failedVectors, remoteMates, flags);

      totalStats.nRead         += regionStats.nRead;
      totalStats.nWritten      += regionStats.nWritten;
      totalStats.nReversed     += regionStats.nReversed;
      totalStats.nOverEnds     += regionStats.nOverEnds;
      totalStats.nRemoteMate   += regionStats.nRemoteMate;
      totalStats.nUnmappedMate += regionStats.nUnmappedMate;
    }
  }

  if (verbosity > 0) {
    printf(" Note in stats below, the numbers recorded represent mappings attempted.\n"
           " For reads which map across multiple destination blocks, the read will be counted multiple times\n");
    printf(" Total reads read in mapped regions           %lld\n", totalStats.nRead);
    printf(" Total reads written in mapped regions        %lld\n", totalStats.nWritten);
    printf(" Total reads where orientation reversed       %lld\n", totalStats.nReversed);
    printf(" Total reads with remotely located mates      %lld\n", totalStats.nRemoteMate);
    printf(" Total reads extending beyond ends of region  %lld (not written)\n", totalStats.nOverEnds);
    printf(" Total reads with unmapped mates              %lld (not written)\n", totalStats.nUnmappedMate);
  }

  samclose(out);
  bam_index_destroy(idx);
  samclose(in);

  if (flags & M_VERIFY) {
    out = samopen(outFName, "rb", 0);
    if (out == 0) {
      fprintf(stderr, "Fail to open BAM file %s\n", outFName);
      return 1;
    }

    verifyBam(out);

    samclose(out);
  }

  if (verbosity > 0) printf("Done\n");
  return 0;
}


// Currently just checks that output is in sorted order
int verifyBam(samfile_t *sam) {
  bam1_t *bs[2];
  int firstind  = 0;
  int secondind = 1;
  int i;

  if (verbosity > 0) printf("Verifying\n");

  for (i=0;i<2;i++) bs[i] = bam_init1();
  
  if (bam_read1(sam->x.bam, bs[0]) <=0 || bs[0]->core.tid < 0) {
    fprintf(stderr,"Failed reading first bam entry from output file in verifyBam\n");
    return 0;
  }

  while (bam_read1(sam->x.bam, bs[secondind]) > 0 && bs[secondind]->core.tid >= 0) {
/*
    if (bs[firstind]->core.tid > bs[secondind]->core.tid || 
        bs[firstind]->core.pos > bs[secondind]->core.pos || 
*/

// Note this is the way the comparison is done in bam_sort.c which is why I do it this way
    if (((uint64_t)bs[firstind]->core.tid<<32|(bs[firstind]->core.pos+1)) > ((uint64_t)bs[secondind]->core.tid<<32|(bs[secondind]->core.pos+1))) {
      fprintf(stderr, "bam entries not in sorted order in output file in verifyBam\n");
      printBam(stderr, bs[firstind],  sam->header);
      printBam(stderr, bs[secondind], sam->header);
        
      return 0;
    }
    firstind  ^= 1;
    secondind ^= 1;
  }

  for (i=0;i<2;i++) bam_destroy1(bs[i]);
  return 1;
}

void printMapping(FILE *fp, Mapping *mapping) {
  fprintf(fp, "%s\t%ld\t%ld\t%s\t%ld\t%ld\t%d\n",Slice_getChrName(mapping->destSlice),
                                        Slice_getChrStart(mapping->destSlice),
                                        Slice_getChrEnd(mapping->destSlice),
                                        Slice_getChrName(mapping->sourceSlice),
                                        Slice_getChrStart(mapping->sourceSlice),
                                        Slice_getChrEnd(mapping->sourceSlice),
                                        mapping->ori);
}

int getPairedMappingFailLists(samfile_t *in, Vector **mappingVectors, Vector ***failedVectorsP, Vector ***remoteMatesP) {
  int i;

  int32_t  curtid = -1;
  Vector  *mappings;
  Mapping *curMapping;
  int      mappingInd = 0;

  long long nRead         = 0;
  long long nRead1        = 0;
  long long nRead2        = 0;
  long long nRead1Read2   = 0;
  long long nNoRead       = 0;
  long long nUnmapped     = 0;
  long long nMateUnmapped = 0; // Note not related to nUnmapped - number where flag indicates mate unmapped in original file
  long long nProperPair   = 0; 
  long long nNoOverlap    = 0;
  long long nFUn          = 0;

  bam1_t *b = bam_init1();

  Vector **unmapped = calloc(in->header->n_targets, sizeof(Vector *));

  Vector **remoteMates = calloc(in->header->n_targets, sizeof(Vector *));

// don't know why I need to check for b->core.tid >= 0 but I do otherwise I get a bam entry with -1 tid
  while (bam_read1(in->x.bam, b) > 0 && b->core.tid >= 0) {
    nRead++;

    if ((b->core.flag & (BAM_FREAD1 | BAM_FREAD2)) == (BAM_FREAD1 | BAM_FREAD2)) {
      nRead1Read2++;
    } else if (b->core.flag & BAM_FREAD1) {
      nRead1++;
    } else if (b->core.flag & BAM_FREAD2) {
      nRead2++;
    } else {
      nNoRead++;
    }

    if (b->core.tid != curtid) {
      curtid = b->core.tid;
//      printf("curtid = %d\n",curtid);
      mappings = mappingVectors[curtid];

      unmapped[curtid] = Vector_new();
      remoteMates[curtid] = Vector_new();

      mappingInd=0;
      if (Vector_getNumElement(mappings)) {
        curMapping = Vector_getElementAt(mappings,mappingInd);
      } else {
        curMapping = NULL;
      }
      if (verbosity > 0) printf("Finding failed pair and remote mate mappings in '%s' (tid = %d)\n",in->header->target_name[curtid],curtid);
      //printf("curMapping = %d\n",curMapping);

      //if (curtid > 1) break;
    }
    
    int end;

    end = bam_calend(&b->core, bam1_cigar(b));

    // Move to a point where the slice ends after the start of the feature (pos)
    // So now either:
    //    Pos starts before the slice and ends before it  - unmapped
    //    Pos starts before the slice and ends in the slice - unmapped
    //    Pos starts before the slice and ends after the slice - unmapped
    //    Pos and end are both in the slice - MAPPED
    //    Pos starts in the slice and ends after it - unmapped
    while (curMapping && b->core.pos > Slice_getChrEnd(curMapping->sourceSlice)-1) {
      if (Vector_getNumElement(mappings) > mappingInd) {
        curMapping = Vector_getElementAt(mappings,mappingInd++);
        if (verbosity > 1) printMapping(stdout, curMapping);
      } else {
        curMapping = NULL;
      }
    }

    if (b->core.flag & BAM_FUNMAP) {
      nFUn++;
    }


    if (b->core.flag & BAM_FMUNMAP) {
      nMateUnmapped++;
    } else {
      // Special case 1 - read lies after end of last assembly mapping block on current sequence, so ignored in later destSlice based fetching
      if (!curMapping) {
        nUnmapped++;
        nNoOverlap++;
      //} else if (b->core.pos >= Slice_getChrStart(curMapping->sourceSlice)-1 &&
      //           end < Slice_getChrEnd(curMapping->sourceSlice)) {
      // Special case 2 - end lies before start of slice, means feature lies completely between mapped regions, so ignored in later destSlice based fetching
      } else if (end < Slice_getChrStart(curMapping->sourceSlice)) {
        nUnmapped++;
        nNoOverlap++;
      // Read which hangs off the end of a mapping block - save for later checks
      } else if (b->core.pos < Slice_getChrStart(curMapping->sourceSlice)-1 || 
                 end > Slice_getChrEnd(curMapping->sourceSlice)) {
        nUnmapped++;
        Vector_addElement(unmapped[curtid],bam_dup1(b));
      }
    }

// I need to store both of pair, because don't know order in which mates will be accessed when doing things by dest slice in main mapping
// Go through each read
//   When start a new seq region, make an array for remote reads 
// Map mate
//   If remote 
//     If on a later seqid, later in this sequence
//       Store read in vector for each seq_region
//     else  // so its the second occuring read of the pair (in seqid, pos terms)
//       Store read in vector for each seq_region
//     endif
//   endif
// At end, sort each vector on pos and name - will be sorted by pos already because came to them in order, so really just subsorting by name

    // Only store reads which map cleanly, so we don't try later to use reads which don't map as mates
    if (curMapping && b->core.pos >= Slice_getChrStart(curMapping->sourceSlice)-1 && 
                 end <= Slice_getChrEnd(curMapping->sourceSlice)) {
      if (b->core.mtid >= 0 && curMapping) {
  
  //   If remote mate
        if (!(b->core.mpos < Slice_getChrEnd(curMapping->sourceSlice) && 
              b->core.mpos >= Slice_getChrStart(curMapping->sourceSlice)-1 && 
              b->core.tid == b->core.mtid)) {
          Mapping *containingMapping;
  
          if ((mapRemoteLocation(mappingVectors, b->core.mtid, b->core.mpos+1, &containingMapping) - 1) < 0) {
  //        Mate doesn't map at all - no need to store
          } else {
  
  // Mate does map (maybe only partially), so need to store
  
  //  Actually now no difference in whats done for each mate so may scrap test - leaving for now
  //      If on a later seqid or later in this sequence
            if (b->core.mtid > b->core.tid || b->core.mpos > b->core.pos) {
  //       Store read in remoteMates vector for curtid
              Vector_addElement(remoteMates[curtid], bam_dup1(b));
              //printBam(stdout, b,in->header);
            } else {
  //     else  // so its the second occuring read of the pair (in seqid, pos terms)
  //       Store read in remoteMates vector for curtid
              Vector_addElement(remoteMates[curtid], bam_dup1(b));
              //printBam(stdout, b,in->header);
            }
          }
        }
      }
    }
  }


  printf("n target = %d\n",in->header->n_targets);
  fflush(stdout);
  for (i=0;i<in->header->n_targets;i++) {
    if (!unmapped[i]) {
      printf("NO VECTOR FOR tid %d (%s)\n",i,in->header->target_name[i]);
    } else {
      int j;
      printf("Number of partially mapped for tid %d (%s) is %d\n",i,in->header->target_name[i],Vector_getNumElement(unmapped[i]));
      if (verbosity > 3) {
        for (j=0;j<Vector_getNumElement(unmapped[i]);j++) {
          bam1_t *b = (bam1_t *) Vector_getElementAt(unmapped[i],j);
          printBam(stdout, b, in->header);
        }
      }
    }
  }

  long long nRemoteMate = 0;
  for (i=0;i<in->header->n_targets;i++) {
    if (remoteMates[i]) {
      Vector_sort(remoteMates[i], bamPosNameCompFunc);
      nRemoteMate += Vector_getNumElement(remoteMates[i]);
    }
  }

  for (i=0;i<in->header->n_targets;i++) {
    if (unmapped[i]) {
      Vector_sort(unmapped[i], bamPosNameCompFunc);
    }
  }

  printf("n target = %d\n",in->header->n_targets);
  fflush(stdout);
  for (i=0;i<in->header->n_targets;i++) {
    if (!remoteMates[i]) {
      printf("NO VECTOR FOR tid %d (%s)\n",i,in->header->target_name[i]);
    } else {
      int j;
      printf("Number of remoteMates for tid %d (%s) is %d\n",i,in->header->target_name[i],Vector_getNumElement(remoteMates[i]));
      if (verbosity > 3) {
        for (j=0;j<Vector_getNumElement(remoteMates[i]);j++) {
          bam1_t *b = (bam1_t *) Vector_getElementAt(remoteMates[i],j);
          printBam(stdout, b, in->header);
        }
      }
    }
  }


  if (verbosity > 0) {
    printf("Total number of reads in input file =         %lld\n",nRead);
    printf("Total READ1 reads in input file =             %lld\n",nRead1);
    printf("Total READ2 reads in input file =             %lld\n",nRead2);
    printf("Total (READ1 & READ2) reads in input file =   %lld\n",nRead1Read2);
    printf("Total !(READ1 | READ2) reads in input file =  %lld\n",nNoRead);
    printf("Total unmapped (FUNMAP) in input file =       %lld\n",nFUn);
    printf("Total mate unmapped (FMUNMAP) in input file = %lld\n\n",nMateUnmapped);
    printf("Total unmapped in assembly mapping =          %lld\n",nUnmapped);
    printf("Total with no overlap with map regions =      %lld\n",nNoOverlap);
    printf("Total number of remote mates =                %lld\n",nRemoteMate);
  }

  *failedVectorsP = unmapped;
  *remoteMatesP   = remoteMates;
  return 1;
}

int bamPosNameCompFunc(const void *one, const void *two) {
/*
  bam1_t **b1p = (bam1_t **)one;
  bam1_t **b2p = (bam1_t **)two;

  bam1_t *b1 = *b1p;
  bam1_t *b2 = *b2p;
*/

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

/*
 Program usage message
*/
void Bammap_usage() {
  printf("bammap \n"
         "  -i --in_file     Input BAM file to map from (string)\n"
         "  -o --out_file    Output BAM file to write (string)\n"
         "  -U --ucsc_naming Input BAM file has 'chr' prefix on ALL seq region names (flag)\n"
         //"  -q --quick       Q&D mapping (flag)\n"
         "  -h --host        Database host name for db containing mapping (string)\n"
         "  -n --name        Database name for db containing mapping (string)\n"
         "  -u --user        Database user (string)\n"
         "  -p --password    Database password (string)\n"
         "  -P --port        Database port (int)\n"
         "  -s --source_ass  Assembly name to map from (string)\n"
         "  -d --dest_ass    Assembly name to map to (string)\n"
         "  -v --verbosity   Verbosity level (int)\n"
         "  -V --verify      Verify output BAM (currently just checks is in sorted order)\n"
         "\n"
         "Notes:\n"
         "  -U will cause 'chr' to be prepended to all source seq_region names, except the special case MT which is changed to chrM.\n"
         "  -v Default verbosity level is 1. You can make it quieter by setting this to 0, or noisier by setting it > 1.\n"
         //"  -q speeds up the process by not pre screening for reads which only partially mapped.\n" 
         //"     This means it can't know when a read has a mate which only partially maps, so can't do the extra\n"
         //"     fixup work on flag, mpos and mtid for these.\n"
         );
  exit(1);
}

/* 
  Take header from input BAM file, replace the seq regions, and then write it out
*/
samfile_t *writeBamHeader(char *inFName, char *outFName, Vector *destinationSlices) {
  bam_header_t *destHeader;
  char          out_mode[5];
  samfile_t *   in;
  samfile_t *   out;
  char          line[1024];
  char *        buff = NULL;
  char **       tokens = NULL;
  int           ntoken;

#ifdef _PBGZF_USE
  bam_set_num_threads_per(1);
#endif
  in = samopen(inFName, "rb", 0);
  if (in == 0) {
    fprintf(stderr, "Fail to open BAM file %s\n", inFName);
    return NULL;
  }

  // For debugging - these lines print input BAM file header to stdout
  if (verbosity > 2) {
    printf("Input BAM file header:\n");
    strcpy(out_mode, "wh");
    samopen("-",out_mode,in->header);
  }

  /* Create the @SQ lines for the destination set of sequences */
  int i;
  for (i=0; i<Vector_getNumElement(destinationSlices); i++) {
    Slice *slice = Vector_getElementAt(destinationSlices,i);
    sprintf(line,"@SQ\tSN:%s\tLN:%ld\n", Slice_getChrName(slice), 
                                        Slice_getChrEnd(slice));
    if (buff) {
      buff = StrUtil_appendString(buff,line);
    } else {
      buff = StrUtil_copyString(&buff,line,0);
    }
  }

  destHeader = bam_header_init();

  /* add non @SQ header lines from source header into destination header */
  //fprintf(stderr,"Tokenize\n");
  StrUtil_tokenizeByDelim(&tokens, &ntoken, in->header->text, "\n");
  //fprintf(stderr,"Token loop\n");
  //fprintf(stderr,"ntoken = %d\n",ntoken);
  for (i=0;i<ntoken;i++) {
    //fprintf(stderr,"token = %s\n",tokens[i]);
    if (strncmp(tokens[i],"@SQ",3)) {
      /* HD must come first */
      if (!strncmp(tokens[i],"@HD",3)) {
        char *tmpBuff;
        StrUtil_copyString(&tmpBuff,tokens[i],0);
        if (!strstr(tokens[i],"SO:coordinate")) {
          tmpBuff = StrUtil_appendString(tmpBuff,"\tSO:coordinate");
        }
        tmpBuff = StrUtil_appendString(tmpBuff,"\n");
        tmpBuff = StrUtil_appendString(tmpBuff,buff);
        free(buff);
        buff = tmpBuff;
      } else {
        //fprintf(stderr,"Here with token = %s buff = %d\n",tokens[i], buff);
        buff = StrUtil_appendString(buff,tokens[i]);
        buff = StrUtil_appendString(buff,"\n");
      }
    }
  }

  //fprintf(stderr,"Before header fix\n");
  // Wasn't an HD line so add one to flag that the file is sorted
  if (strncmp(buff,"@HD",strlen("@HD"))) {
    char *tmpBuff;
    //fprintf(stderr,"Header fix\n");
    StrUtil_copyString(&tmpBuff,"@HD\tVN:1.0\tSO:coordinate\n",0);
    tmpBuff = StrUtil_appendString(tmpBuff,buff);
    free(buff);
    buff = tmpBuff;
  }

  destHeader->l_text = strlen(buff);
  destHeader->text   = buff;
  //fprintf(stderr, "header->text = %s\n",buff);
  sam_header_parse(destHeader);

  strcpy(out_mode, "wbh");

#ifdef _PBGZF_USE
  bam_set_num_threads_per(9);
#endif
  out = samopen(outFName,out_mode,destHeader);

/* Need this hash initialised for looking up tids */
  bam_init_header_hash(out->header);

  if (verbosity > 2) {
    printf("Output BAM file header:\n");
    strcpy(out_mode, "wh");
    samopen("-",out_mode,out->header);
  }

  samclose(in);

  return out;
}

void clearPairing(bam1_t *b, int end) {
  b->core.mpos = -1;
  b->core.mtid = -1;

  b->core.isize = end - b->core.pos;

  b->core.flag |= BAM_FMUNMAP;
  b->core.flag &= ~(BAM_FPROPER_PAIR);
}

int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 9, 14, 1, 6, 5, 13, 3, 11, 7, 15 };

int mapBam(char *fName, samfile_t *out, Mapping *mapping, ReadMapStats *regionStats, 
           samfile_t *in, bam_index_t *idx, Vector **mappingVectors, Vector **failedVectors, Vector **remoteMates, int flags) {
  int  ref;
  int  begRange;
  int  endRange;
  char region[1024];


  sprintf(region,"%s:%ld-%ld", Slice_getChrName(mapping->sourceSlice), 
                             Slice_getChrStart(mapping->sourceSlice), 
                             Slice_getChrEnd(mapping->sourceSlice));
  bam_parse_region(in->header, region, &ref, &begRange, &endRange);
  if (ref < 0) {
    fprintf(stderr, "Invalid region %s %ld %ld\n", Slice_getChrName(mapping->sourceSlice), 
                                                 Slice_getChrStart(mapping->sourceSlice), 
                                                 Slice_getChrEnd(mapping->sourceSlice));
    return 1;
  }

  bam_iter_t iter = bam_iter_query(idx, ref, begRange, endRange);
  bam1_t *b = bam_init1();

  int32_t newtid = bam_get_tid(out->header, Slice_getChrName(mapping->destSlice));

  if (mapping->ori == 1) {
    mapBam_forward(mapping, b, in, out, mappingVectors, failedVectors, remoteMates, begRange, endRange, regionStats, newtid, iter);
  } else if (mapping->ori == -1) {
    mapBam_reverse(mapping, b, in, out, mappingVectors, failedVectors, remoteMates, begRange, endRange, regionStats, newtid, iter);
  } else {
    fprintf(stderr,"Error: Unknown mapping orientation %d\n", mapping->ori);
  }

  bam_iter_destroy(iter);
  bam_destroy1(b);
}


int mapBam_reverse(Mapping *mapping, bam1_t *b, samfile_t *in, samfile_t *out, Vector **mappingVectors, Vector **failedVectors, Vector **remoteMates, 
                   int begRange, int endRange, ReadMapStats *regionStats, int newtid, bam_iter_t iter) {
  Vector *reverseCache = Vector_new();
  Vector *reverseFinal = Vector_new();
  int i;

  while (bam_iter_read(in->x.bam, iter, b) >= 0) {
    Vector_addElement(reverseCache,bam_dup1(b));
  }
 
  Vector_sort(reverseCache, bamPosNameCompFunc);
    
  if (verbosity > 2) { printf("reverseCache with %d elements\n",Vector_getNumElement(reverseCache)); }

  for (i=0;i < Vector_getNumElement(reverseCache); i++) { 
    int end;

    // Hacky - for reverse need to make a copy so we don't modify the one in the reverseCache array which is used for mate finding
    bam_copy1(b, Vector_getElementAt(reverseCache, i));
    // For the copy, we clear the MY_FUSEDFLAG which can only have been added by the mate finding code to entry in reverseCache
    // We DON'T want to write that flag
    b->core.flag &= ~(MY_FUSEDFLAG);

    regionStats->nRead++;

    end = bam_calend(&b->core, bam1_cigar(b));

    // There is a special case for reads which have zero length and start at begRange (so end at begRange ie. before the first base we're interested in).
    // That is the reason for the || end == begRange test
    if (end > endRange || b->core.pos < begRange || end == begRange) {
      regionStats->nOverEnds++;
      continue;
    }
    
    // pos should be left most position on reference so is either pos or end depending on mapping orientation
    int newpos = calcNewPos(mapping, b, end);

    if (b->core.mtid >= 0) {
      mapMateLocation(mapping, b, end, in, out, mappingVectors, failedVectors, remoteMates, reverseCache, regionStats, begRange, endRange, newtid, newpos);
    }

    b->core.tid = newtid;
    b->core.pos = newpos;

    // Extra stuff needed for a reverse mapping
    // toggle (XOR) rev com flag if reverse orientation mapping
    // Also have to revcom the sequence
    // Also have to reverse the quality
    // Also have to reverse the cigar
    b->core.flag ^= BAM_FREVERSE;

    regionStats->nReversed++;

    uint8_t *seq;
    uint8_t *qual;
    seq = bam1_seq(b);
    qual = bam1_qual(b);

    int i;
    uint8_t *revseq  = calloc(1, (b->core.l_qseq+1)/2 + b->core.l_qseq);
    uint8_t *revqual = revseq + (b->core.l_qseq+1)/2;
    
    // make rev com seq and reverse qual at same time
    for (i=0;i<b->core.l_qseq;i++) {
      uint8_t ib = seq_comp_table[bam1_seqi(seq, b->core.l_qseq-i-1)];
      revseq[i/2] |= ib << 4*(1-i%2);
      revqual[i] = qual[b->core.l_qseq-i-1];
    }

    memcpy(seq, revseq, (b->core.l_qseq+1)/2 + b->core.l_qseq);
    free(revseq);
    
    int *revcig = calloc(1,b->core.n_cigar*4);
    int *cig = bam1_cigar(b);
    
    for (i=0;i<b->core.n_cigar;i++) {
      revcig[i] = cig[b->core.n_cigar-1-i];
    }
    memcpy(cig,revcig,b->core.n_cigar * 4);

    free(revcig);

    // recalculate bin, not sure if need this for all coords, or just - ori ones
    // Note: Here we are working with the mapped position (pos has been changed)
    end = bam_calend(&b->core, bam1_cigar(b));
    b->core.bin = bam_reg2bin(b->core.pos,end);
  
    // For reverse rather than writing straight away, store the transformed entries, so can sort into correct order before writing
    Vector_addElement(reverseFinal, bam_dup1(b));
  }

  Vector_sort(reverseFinal, bamPosCompFunc);
  for (i=0;i<Vector_getNumElement(reverseFinal);i++) {
    regionStats->nWritten++;

    if (!bam_write1(out->x.bam, Vector_getElementAt(reverseFinal,i))) {
      fprintf(stderr, "Failed writing bam entry\n");
    }
  }

  // Need to free the cache
  // Note bam_destroy1 isn't a function, its a macro, so can't be used as a function pointer for Vector free function
  for (i=0;i<Vector_getNumElement(reverseCache);i++) {
    bam_destroy1((bam1_t *)Vector_getElementAt(reverseCache,i));
  }
  for (i=0;i<Vector_getNumElement(reverseFinal);i++) {
    bam_destroy1((bam1_t *)Vector_getElementAt(reverseFinal,i));
  }
    
  Vector_free(reverseCache);
  Vector_free(reverseFinal);

  return 0;
}

int mapBam_forward(Mapping *mapping, bam1_t *b, samfile_t *in, samfile_t *out, Vector **mappingVectors, Vector **failedVectors, Vector **remoteMates, 
                   int begRange, int endRange, ReadMapStats *regionStats, int newtid, bam_iter_t iter) {
  while (bam_iter_read(in->x.bam, iter, b) >= 0) {
    int end;

  //fprintf(stderr, "HERE3\n");
    regionStats->nRead++;

    end = bam_calend(&b->core, bam1_cigar(b));

    // There is a special case for reads which have zero length and start at begRange (so end at begRange ie. before the first base we're interested in).
    // That is the reason for the || end == begRange test
    if (end > endRange || b->core.pos < begRange || end == begRange) {
      regionStats->nOverEnds++;
      continue;
    }
    
    // pos should be left most position on reference so is either pos or end depending on mapping orientation
    int newpos = calcNewPos(mapping, b, end);

    if (b->core.mtid >= 0) {
      mapMateLocation(mapping, b, end, in, out, mappingVectors, failedVectors, remoteMates, NULL, regionStats, begRange, endRange, newtid, newpos);
    }

    b->core.tid = newtid;
    b->core.pos = newpos;

    // recalculate bin, not sure if need this for all coords, or just - ori ones
    // Note: Here we are working with the mapped position (pos has been changed)
    end = bam_calend(&b->core, bam1_cigar(b));
    b->core.bin = bam_reg2bin(b->core.pos,end);

    regionStats->nWritten++;

    if (!bam_write1(out->x.bam, b)) {
      fprintf(stderr, "Failed writing bam entry\n");
    }
  }

  return 0;
}


int mapMateLocation(Mapping *mapping, bam1_t *b, int end, samfile_t *in, samfile_t *out, Vector **mappingVectors, 
                    Vector **failedVectors, Vector **remoteMates, Vector *reverseCache, 
                    ReadMapStats *regionStats, int begRange, int endRange, int newtid, int newpos) {
// Mate handling is by far the most complicated part of this process. In summary:
//    Check if the mate only partially maps (is in failedVector). 
//       If not then clear mate information for this read and recalculate isize
//    else it does map cleanly or doesn't map at all!
//       If its remote
//         If it doesn't map at all
//           Clear mate information for this read and recalculate isize
//         else (so it maps)
//           Check if we can not find its mate in remoteMates
//             Clear mate information for this read and recalculate isize
//           else (so can find mate)
//             Use mate to calculate new isize etc
//       else (so its local)
//         Is it forward strand
//           If it doesn't map at all ?? Not sure this can happen
//             Clear mate information for this read and recalculate isize
//           else (so it maps)
//             Nothing extra to do
//         Is it reverse strand
//           If it doesn't map at all ?? Not sure this can happen
//             Clear mate information for this read and recalculate isize
//           else (so it maps)
//             Check if we can not find its mate in reverseCache
//               Clear mate information for this read and recalculate isize
//             else (so can find mate)
//               Use mate to calculate new isize etc
//           
  bam1_t *fb;
  if ((fb = mateFoundInVectors(b, failedVectors))) {
    clearPairing(b, end);
    // Eliminate from search of failed array by setting one of the unused BAM flag bits
    // Reason for doing this is that in complex mapping cases there can be multiple mappings at same
    // location so don't want to eliminate all of them if one mate doesn't map but others do
    fb->core.flag |= MY_FUSEDFLAG;
  } else {
    // If mate lies outside current mapping block
    if (!(b->core.mpos < endRange && b->core.mpos >= begRange && b->core.tid == b->core.mtid)) {
      Mapping *containingMapping;
      int newmpos;
  
      regionStats->nRemoteMate++;
      if ((newmpos = mapRemoteLocation(mappingVectors, b->core.mtid, b->core.mpos+1, &containingMapping) - 1) < 0) {
        regionStats->nUnmappedMate++;
  
        clearPairing(b, end);
      } else {
        //printf("containing slice = %s %d %d\n",Slice_getChrName(containingMapping->destSlice),
        //                                       Slice_getChrStart(containingMapping->destSlice),
        //                                       Slice_getChrEnd(containingMapping->destSlice));
        int32_t newmtid = bam_get_tid(out->header, Slice_getChrName(containingMapping->destSlice));

        //printf("Searching for mate of:\n");
        //printBam(stdout,b,in->header);
        bam1_t *mate = getMateFromRemoteMates(b, remoteMates);

        if (!mate) {
          // Shouldn't happen, all remote mates should be in remoteMates
          fprintf(stderr,"Error: Missing remote mate for %s. Shouldn't happen\n", bam1_qname(b));
          printBam(stderr,b,in->header);

          clearPairing(b, end);
        } else {
          // use mate to set everything correctly, isize etc!
          // If mapping block containing mate is - ori then need to switch 'strand' of mpos
          //    'end' position of mate read will now be mpos (because of - ori mapping)

          // Actually some depends on RELATIVE orientation of two mapping blocks
          //    If both mapping and containingMapping are +
          //    If both mapping and containingMapping are -
          //    If mapping is - and containingMapping is +
          //    If mapping is + and containingMapping is -

          // newpos is already corrected for orientation of mapping

          // Note: We can't freely use mapLocation on mate end position because although we know that
          //       mate isn't in failed mates vector, if the -q flag is used then that vector is not filled
          //       So first thing - check that mate end also maps in containingMapping
          int tmend = bam_calend(&mate->core, bam1_cigar(mate));

          if (tmend > Slice_getChrEnd(containingMapping->sourceSlice)) {
            fprintf(stderr,"Error: Remote mate doesn't cleanly map for %s. Shouldn't happen\n", bam1_qname(b));
            fprintf(stderr,"  containing slice = %s %ld %ld\n",
                    Slice_getChrName(containingMapping->sourceSlice),
                    Slice_getChrStart(containingMapping->sourceSlice),
                    Slice_getChrEnd(containingMapping->sourceSlice)
                   );
            fprintf(stderr,"  b   : ");
            printBam(stderr, b,in->header);
            fprintf(stderr,"  mate: ");
            printBam(stderr, mate,in->header);

            clearPairing(b, end);
          } else {
            // First fix newmpos which may in fact be new end position
            if (containingMapping->ori == -1) {
              newmpos  = calcNewPos(containingMapping, mate, tmend);
            }

            // If mates both map to same sequence, calculate new isize
            if (newmtid == newtid) {
              // To simplify the logic I'm going to calculate all the mapped positions and flags given the orientation of each mapping block
              // and then recalculate isize 
              int bpos, bend, matepos, mateend;
              int bflag = 0, mateflag = 0;
              uint32_t b5, mate5;
         
              matepos = newmpos;  // newmpos is already corrected for orientation
              mateend = calcNewEnd(containingMapping, mate, tmend);
              if (containingMapping->ori == 1) {
                mateflag = mate->core.flag;
              } else {
                mateflag = mate->core.flag^BAM_FREVERSE;
              }
  
              bpos = newpos; // newpos is already corrected for orientation
              bend = calcNewEnd(mapping, b, end);
              if (mapping->ori == 1) {
                bflag = b->core.flag;
              } else {
                bflag = b->core.flag^BAM_FREVERSE;
              }

              b5    = (bflag&BAM_FREVERSE) ? bend : bpos;
              mate5 = (mateflag&BAM_FREVERSE) ? mateend : matepos;
  
              //printf("Old isize = %d new isize = %d\n", b->core.isize, (mate5-b5));
              b->core.isize = mate5 - b5;
            } else {
              b->core.isize = 0;
            }
  
            b->core.mpos = newmpos;
            b->core.mtid = newmtid;

            if (containingMapping->ori == -1) {
              b->core.flag ^= BAM_FMREVERSE;
            }

            // If mates are on different orientation mapping blocks, proper pairing rule likely to be broken so clear proper pair flag
            if (containingMapping->ori != mapping->ori) {
              b->core.flag &= ~(BAM_FPROPER_PAIR);
            }

            // mark mate as used
            mate->core.flag |= MY_FUSEDFLAG;
          }
        }
      }
    } else if (mapping->ori == 1) { // For local, forward strand mapping
      b->core.mtid = newtid;
  
      if ((b->core.mpos = mapLocation(mapping, b->core.mpos+1) - 1) < 0) {
        regionStats->nUnmappedMate++;

        clearPairing(b, end);
      }
    } else { // local reverse mapping
             // have local cache of reads (reverseCache) so look up mate in there
             // Change mpos on b to end of mate
             // Reverse FMREVERSE
             // Mark mate as used
      b->core.mtid = newtid;
  
      bam1_t *mate = findMateInVector(b, reverseCache);

      if (!mate) {
        // Shouldn't happen, all local reverse mates should be in rreverseMates
        fprintf(stderr,"Error: Missing local reverse mate for %s. Shouldn't happen\n", bam1_qname(b));
        printBam(stderr,b,in->header);
  
        clearPairing(b, end);
      } else {
        int mend = bam_calend(&mate->core, bam1_cigar(mate));
        if (mend > Slice_getChrEnd(mapping->sourceSlice)) {
          fprintf(stderr,"Error: Reverse local mate doesn't cleanly map for %s. Shouldn't happen\n", bam1_qname(b));
          fprintf(stderr,"  containing slice = %s %ld %ld\n",
                  Slice_getChrName(mapping->sourceSlice),
                  Slice_getChrStart(mapping->sourceSlice),
                  Slice_getChrEnd(mapping->sourceSlice)
                 );
          fprintf(stderr,"  b   : ");
          printBam(stderr, b,in->header);
          fprintf(stderr,"  mate: ");
          printBam(stderr, mate,in->header);

          clearPairing(b, end);
        } else {
          b->core.mpos = calcNewPos(mapping, mate, bam_calend(&mate->core, bam1_cigar(mate)));
          b->core.flag ^= BAM_FMREVERSE;
    
          mate->core.flag |= MY_FUSEDFLAG;
        }
      }
    }
  }
}

int calcNewEnd(Mapping *mapping, bam1_t *b, int sourceEnd) {
  int newpos;

  if (mapping->ori == 1) {
    newpos = mapLocation(mapping, sourceEnd) - 1;
  } else {
    newpos = mapLocation(mapping, b->core.pos+1) - 1;
  }

  return newpos;
}

int calcNewPos(Mapping *mapping, bam1_t *b, int sourceEnd) {
  int newpos;

  if (mapping->ori == 1) {
    newpos = mapLocation(mapping, b->core.pos+1) - 1;
  } else {
    newpos = mapLocation(mapping, sourceEnd) - 1;
  }

  return newpos;
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

int mapRemoteLocation(Vector **mappingVectors, int seqid, int pos, Mapping **containingMappingP) {
  Mapping *mapping = NULL;
  Vector  *mapVec = mappingVectors[seqid];
  int i;
  int imin = 0;
  int imax = Vector_getNumElement(mapVec)-1;

  *containingMappingP = NULL;

  /* Binary search to find the mapping containing the location */

  while (imax >= imin) {
    int imid = (imax+imin) / 2;
 
    Mapping *m = Vector_getElementAt(mapVec,imid);
//    printf("imid = %d imin = %d imax = %d pos = %d m start = %d m end = %d\n",imid,imin,imax,pos,Slice_getChrStart(m->sourceSlice),Slice_getChrEnd(m->sourceSlice));

    if (pos > Slice_getChrEnd(m->sourceSlice)) {
      imin = imid + 1;
    } else if (pos < Slice_getChrStart(m->sourceSlice)) {
      imax = imid - 1;
    } else {
      // key found at index imid
      mapping = m;
      break;
    }
  }

  if (!mapping) {
    // fprintf(stderr,"Mate location lies outside any mapped region pos = %d\n", pos);
    return -1;
  }

  *containingMappingP = mapping;
  return mapLocation(mapping, pos);
}

inline int mapLocation(Mapping *mapping, int pos) {
  int fromStart = Slice_getChrStart(mapping->sourceSlice);
  int fromEnd   = Slice_getChrEnd(mapping->sourceSlice);
  int toStart   = Slice_getChrStart(mapping->destSlice);
  int toEnd     = Slice_getChrEnd(mapping->destSlice);

  if (pos < fromStart || pos > fromEnd) {
    fprintf(stderr,"Error: tried to map position out of range pos = %d range = %s %d %d\n", pos, Slice_getChrName(mapping->sourceSlice), fromStart, fromEnd);
    return -1;
  }

  if (mapping->ori == 1) {
    pos += (toStart-fromStart);
  } else {
    int offset = pos-fromStart;
    pos = toEnd-offset;
  }

  return pos;
}


Vector **getMappingVectorsBySourceRegion(DBAdaptor *dba, samfile_t *in, char *sourceName, char *destName, int flags) {
  int i;
  int str_offset = 0;
  char *MT = "MT";

  if (flags & M_UCSC_NAMING) str_offset = 3;

  Vector **mappingVectors = calloc(in->header->n_targets, sizeof(Vector *));

  for (i=0;i<in->header->n_targets;i++) {
// For UCSC naming chop off 'chr' prefix
    char *seqName = &(in->header->target_name[i][str_offset]);

// UCSC MT
    if (flags & M_UCSC_NAMING && !strcmp(seqName,"M")) seqName = MT;
      

    // Note reverse mapping direction to key on source
    mappingVectors[i] = getMappings(dba,seqName,destName,sourceName, 1, flags); 
  }

  return mappingVectors;
}

// rev flag is a bit of a hack, to enable fetching by source region, but switching slices so
// it looks like the mapping is the other way round 
Vector *getMappings(DBAdaptor *dba, char *seqName, char *fromAssName, char *toAssName, int rev, int flags) {
  StatementHandle *sth;
  ResultRow *      row;
  char             qStr[1024];
  char *           coordSysName1;
  char *           coordSysName2;
  Vector *         mappingVector = Vector_new();

  int i;
  for (i=0;i<2;i++) {
    if (!i) {
      coordSysName1 = toAssName;
      coordSysName2 = fromAssName;
    } else {
      coordSysName1 = fromAssName;
      coordSysName2 = toAssName;
    }
    sprintf(qStr,
            "select sr1.name, a.asm_start, a.asm_end, cs1.coord_system_id, sr1.length, sr2.name, a.cmp_start, a.cmp_end, cs2.coord_system_id, sr2.length, a.ori"
            " from seq_region sr1, seq_region sr2, assembly a, coord_system cs1, coord_system cs2"
            " where sr1.seq_region_id=a.asm_seq_region_id and"
            "       sr2.seq_region_id=a.cmp_seq_region_id and"
            "       cs1.coord_system_id=sr1.coord_system_id and"
            "       cs2.coord_system_id=sr2.coord_system_id and"
            "       cs1.version='%s' and"
            "       cs2.version='%s' and"
            "       sr%d.name='%s' order by a.%s_start",
            coordSysName1, coordSysName2, i+1,  seqName, (!i ? "asm" : "cmp"));

    if (verbosity > 2) printf("%s\n",qStr);
  
    sth = dba->dbc->prepare(dba->dbc,qStr,strlen(qStr));
  
    sth->execute(sth);
  
    CoordSystemAdaptor *csa = DBAdaptor_getCoordSystemAdaptor(dba);
    SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(dba);
  
    while ((row = sth->fetchRow(sth))) {
      char *destName;
      int   destStart;
      int   destEnd;
      int   destLen;
      char *sourceName;
      int   sourceStart;
      int   sourceEnd;
      int   sourceLen;
      int   ori;
      IDType destCsId;
      IDType sourceCsId;
      Mapping *mapping;
      Slice *dest;
      Slice *source;
  
      if (row->col(row,0))  destName   = row->getStringAt(row,0);
      if (row->col(row,1))  destStart  = row->getIntAt(row,1);
      if (row->col(row,2))  destEnd    = row->getIntAt(row,2);
      if (row->col(row,3))  destCsId   = row->getLongLongAt(row,3);
      if (row->col(row,4))  destLen    = row->getIntAt(row,4);
        
      if (row->col(row,5))  sourceName   = row->getStringAt(row,5);
      if (row->col(row,6))  sourceStart  = row->getIntAt(row,6);
      if (row->col(row,7))  sourceEnd    = row->getIntAt(row,7);
      if (row->col(row,8))  sourceCsId   = row->getLongLongAt(row,8);
      if (row->col(row,9))  sourceLen    = row->getIntAt(row,9);
      if (row->col(row,10))  ori         = row->getIntAt(row,10);

// Maybe should do this after getting mappings - for now do it here
      if (flags & M_UCSC_NAMING) {
        char **chP;
        char  tmp[1024];
        
        if (!rev) {
          if (!i) { chP = &sourceName; }
          else { chP = &destName; }
        } else {
          if (!i) { chP = &destName; }
          else { chP = &sourceName; }
        }

        // Special case for MT -> chrM
        if (!strcmp(*chP,"MT")) (*chP)[1] = '\0';

        strcpy(tmp,"chr");
        *chP = strcat(tmp,*chP);
      }
       
      CoordSystem *destCs = CoordSystemAdaptor_fetchByDbID(csa, destCsId);
      dest    = Slice_new(destName,destStart,destEnd,1,destLen,destCs,sa);
      CoordSystem *sourceCs = CoordSystemAdaptor_fetchByDbID(csa, sourceCsId);
      source  = Slice_new(sourceName,sourceStart,sourceEnd,1,sourceLen, sourceCs,sa);
      mapping = calloc(1,sizeof(Mapping)); 
  
      if (!rev) {
        if (!i) {
          mapping->sourceSlice = source;
          mapping->destSlice   = dest;
          mapping->ori         = ori;
        } else {
          mapping->sourceSlice = dest;
          mapping->destSlice   = source;
          mapping->ori         = ori;
        }
      } else {
        if (!i) {
          mapping->sourceSlice = dest;
          mapping->destSlice   = source;
          mapping->ori         = ori;
        } else {
          mapping->sourceSlice = source;
          mapping->destSlice   = dest;
          mapping->ori         = ori;
        }
      }
      
      Vector_addElement(mappingVector,mapping); 
    }
  
    sth->finish(sth);
  }

  return mappingVector;
}

Vector *getDestinationSlices(DBAdaptor *dba, char *assName) {
  StatementHandle *sth;
  ResultRow *      row;
  char             qStr[1024];
  char *           assemblyType = DBAdaptor_getAssemblyType(dba);

  if (!strcmp(assemblyType,assName)) {
    sprintf(qStr,
            "select distinct(sr.name), sr.length, sr.coord_system_id"
            " from seq_region sr,"
            "      seq_region_attrib sra,"
            "      attrib_type at,"
            "      coord_system cs "
            " where cs.version = '%s' and"
            "       cs.coord_system_id=sr.coord_system_id and"
            "       sra.seq_region_id=sr.seq_region_id and"
            "       sra.attrib_type_id=at.attrib_type_id and"
            "       at.code='toplevel' and"
            "       sr.seq_region_id not in "
            "            (select sr1.seq_region_id from seq_region sr1, seq_region_attrib sra1, attrib_type at1"
            "             where sr1.seq_region_id=sra1.seq_region_id and"
            "                   sra1.attrib_type_id=at1.attrib_type_id and"
            "                   at1.code = 'non_ref') order by length desc",
            assName);
  } else {
    sprintf(qStr,
            "select distinct(sr.name), sr.length, sr.coord_system_id"
            " from seq_region sr,"
            "      coord_system cs "
            " where cs.version = '%s' and"
            "       cs.coord_system_id=sr.coord_system_id order by sr.length desc",
            assName);
  }

  if (verbosity > 2) printf("%s\n",qStr);
  sth = dba->dbc->prepare(dba->dbc,qStr,strlen(qStr));

  sth->execute(sth);

  Vector *toplevelSliceVector = Vector_new();

  CoordSystemAdaptor *csa = DBAdaptor_getCoordSystemAdaptor(dba);
  SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(dba);
  while ((row = sth->fetchRow(sth))) {
    char *name;
    int   length;
    IDType csId;
    Slice *slice;

    if (row->col(row,0))  name   = row->getStringAt(row,0);
    if (row->col(row,1))  length = row->getIntAt(row,1);
    if (row->col(row,2))  csId   = row->getLongLongAt(row,2);

    CoordSystem *cs = CoordSystemAdaptor_fetchByDbID(csa, csId);

     
    slice = Slice_new(name,1,length,1,length,cs,sa);

    Vector_addElement(toplevelSliceVector,slice); 
  }


  sth->finish(sth);

  return toplevelSliceVector;
}
