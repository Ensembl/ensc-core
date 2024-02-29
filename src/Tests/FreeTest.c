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

#include <stdio.h>

#include "SliceAdaptor.h"
#include "DBAdaptor.h"
#include "EnsC.h"
#include "DNAAlignFeature.h"
#include "DNAAlignFeatureAdaptor.h"

#include "BaseRODBTest.h"
#ifdef HAVE_LIBTCMALLOC
#include "gperftools/tcmalloc.h"
#endif

int main(int argc, char *argv[]) {
  DBAdaptor *dba;
  DNAAlignFeatureAdaptor *dafa;
  Slice *slice;
  Vector *features;
  int i;
  int failed = 0;
  int testResult = 0;
  
  initEnsC(argc, argv);

  dba = Test_initROEnsDB();

  slice = Test_getStandardSlice(dba);

  testResult += ok(1, slice!=NULL);

  dafa = DBAdaptor_getDNAAlignFeatureAdaptor(dba);

  testResult += ok(2, dafa!=NULL);

  features =  Slice_getAllDNAAlignFeatures(slice,NULL,NULL,NULL,NULL);

  testResult += ok(3, features!=NULL);
  testResult += ok(4, Vector_getNumElement(features)!=0);

  unsigned long long totNameLen=0;
  unsigned long long totHitNameLen=0;
  unsigned long long totCigarLen=0;
  unsigned long long nFeat = Vector_getNumElement(features);

  for (i=0;i<Vector_getNumElement(features) && !failed;i++) {
    DNAAlignFeature *daf = Vector_getElementAt(features,i);
    Vector *ungapped;
    char *oldCigar = DNAAlignFeature_getCigarString(daf);

    ungapped = DNAAlignFeature_getUngappedFeatures((BaseAlignFeature *)daf);
    if (!ungapped) failed = 1;

    totCigarLen += strlen(DNAAlignFeature_getCigarString(daf));
    totNameLen += strlen(DNAAlignFeature_getSeqName(daf));
    totHitNameLen += strlen(DNAAlignFeature_getHitSeqName(daf));


    BaseAlignFeature_parseFeatures(daf,ungapped); 

    // This is how it was written - not sure about how do do this
    //Object_incRefCount(ungapped);
    //Vector_free(ungapped);
    // I know these will free it but not sure if its best way
    Vector_setFreeFunc(ungapped, Object_freeImpl);
    Vector_free(ungapped);
  }
  fprintf(stderr, "Average hit name len = "IDFMTSTR"\n", totHitNameLen/nFeat);
  fprintf(stderr, "Average name len = "IDFMTSTR"\n", totNameLen/nFeat);
  fprintf(stderr, "Average cigar len = "IDFMTSTR"\n", totCigarLen/nFeat);
  testResult += ok(5, !failed);

#ifdef HAVE_LIBTCMALLOC
  printf("Before calling Vector_free on features\n");
  tc_malloc_stats();
#endif

  Vector_setFreeFunc(features, Object_freeImpl);
  DNAAlignFeatureAdaptor_clearCache(dafa);

#ifdef HAVE_LIBTCMALLOC
  MallocExtension_ReleaseFreeMemory();
  printf("After calling DAFA_cacheClear\n");
  tc_malloc_stats();
#endif

  printf(" I THINK THERE's SOMETHING WRONG IN THE WAY THE ABOVE FREEs AND CLEARS ARE WORKING\n");
  return testResult;
}
