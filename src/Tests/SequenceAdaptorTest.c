/*
 * Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

#include <stdio.h>

#include "SliceAdaptor.h"
#include "SequenceAdaptor.h"
#include "ProjectionSegment.h"
#include "DBAdaptor.h"
#include "EnsC.h"

#include "BaseTest.h"


static char *CHR = "20";
static long START = 30220000;
static long END   = 31200000;
static long STRAND = 1;

void compareComplements(Slice *slice, SequenceAdaptor *seqA);

int nTest = 1;

int main(int argc, char *argv[]) {
  DBAdaptor *dba;
  Slice *slice;
  SliceAdaptor *sliceA;
  SequenceAdaptor *seqA;

  initEnsC(argc, argv);

  //dba = DBAdaptor_new("ensembldb.ensembl.org","anonymous",NULL,"homo_sapiens_core_70_37",5306,NULL);
  //dba = DBAdaptor_new("ens-livemirror.internal.sanger.ac.uk","ensro",NULL,"homo_sapiens_core_70_37",3306,NULL);
  dba = DBAdaptor_new("genebuild2.internal.sanger.ac.uk","ensadmin","ensembl","steve_hs_testdb",3306,NULL);

  sliceA = DBAdaptor_getSliceAdaptor(dba);
  seqA = DBAdaptor_getSequenceAdaptor(dba);

//
// Test fetch_by_Slice_start_end_strand
//
  slice = SliceAdaptor_fetchByRegion(sliceA,"chromosome",CHR,START,END,STRAND,NULL,0);
  compareComplements(slice, seqA);

  slice = SliceAdaptor_fetchByRegion(sliceA, "clone","AL031658.11", POS_UNDEF, POS_UNDEF, STRAND_UNDEF, NULL, 0);
  compareComplements(slice, seqA);

  slice = SliceAdaptor_fetchByRegion(sliceA, "supercontig","NT_028392", POS_UNDEF, POS_UNDEF, STRAND_UNDEF, NULL, 0);
  compareComplements(slice, seqA);

  slice = SliceAdaptor_fetchByRegion(sliceA, "contig", "AL031658.11.1.162976", POS_UNDEF, POS_UNDEF, STRAND_UNDEF, NULL, 0);
  compareComplements(slice, seqA);

  return 0;
}

void compareComplements(Slice *slice, SequenceAdaptor *seqA) {

  char *seq = SequenceAdaptor_fetchBySliceStartEndStrand(seqA, slice,1,POS_UNDEF,1);

  fprintf(stderr,"FORWARD STRAND SLICE SEQ for %s\n", Slice_getName(slice));
  //fprintf(stderr,"%s\n", seq);

  char *invertedSeq = SequenceAdaptor_fetchBySliceStartEndStrand(seqA, slice,1,POS_UNDEF,-1);

  fprintf(stderr,"REVERSE STRAND SLICE SEQ for %s\n", Slice_getName(slice));
  //fprintf(stderr,"%s\n", invertedSeq);

  ok(nTest++, strlen(seq) == Slice_getLength(slice)); //sequence is correct length

  SeqUtil_reverseComplement(seq, Slice_getLength(slice));  //reverse complement seq
  //fprintf(stderr,"%s\n", seq);

  ok(nTest++, !strcmp(seq, invertedSeq)); //revcom same as seq on inverted slice
}

