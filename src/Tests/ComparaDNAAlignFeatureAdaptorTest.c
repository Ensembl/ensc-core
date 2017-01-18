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

#include <stdio.h>

#include "SliceAdaptor.h"
#include "ComparaDBAdaptor.h"
#include "EnsC.h"
#include "ComparaDNAAlignFeatureAdaptor.h"

#include "BaseComparaDBTest.h"

int main(int argc, char *argv[]) {
  ComparaDBAdaptor *cdba;
  ComparaDNAAlignFeatureAdaptor *cdafa;
  Slice *slice = NULL;
  Vector *dnaAligns;
  int i;
  
  initEnsC(argc, argv);

  cdba = Test_initComparaDB();

  slice = Test_getStandardSlice(cdba);

  ok(1, slice!=NULL);

  cdafa = ComparaDBAdaptor_getComparaDNAAlignFeatureAdaptor(cdba);

  ok(2, cdafa!=NULL);

  dnaAligns = ComparaDNAAlignFeatureAdaptor_fetchAllBySlice(cdafa,slice,"mus musculus","NCBIM30","WGA");

  ok(3, dnaAligns!=NULL);
  ok(4, Vector_getNumElement(dnaAligns)!=0);

  for (i=0; i<Vector_getNumElement(dnaAligns); i++) {
    DNAAlignFeature *daf = Vector_getElementAt(dnaAligns,i);
    Slice *slice = (Slice *)DNAAlignFeature_getSlice(daf);
    printf(" %s %ld %ld and %s %d %d\n", DNAAlignFeature_getSeqName((SeqFeature*)daf),
                                       DNAAlignFeature_getStart(daf),
                                       DNAAlignFeature_getEnd(daf),
                                       DNAAlignFeature_getHitSeqName(daf),
                                       DNAAlignFeature_getHitStart(daf),
                                       DNAAlignFeature_getHitEnd(daf));
    printf(" seq = %s\n",Slice_getSubSeq(slice,
                                         DNAAlignFeature_getStart(daf),
                                         DNAAlignFeature_getEnd(daf), 
                                         DNAAlignFeature_getStrand(daf)));
  }

  return 0;
}
