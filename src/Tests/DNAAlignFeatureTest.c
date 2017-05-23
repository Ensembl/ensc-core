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

#include <stdio.h>

#include "SliceAdaptor.h"
#include "DNAAlignFeatureAdaptor.h"
#include "DBAdaptor.h"
#include "EnsC.h"
#include "DNAAlignFeature.h"

#include "BaseRODBTest.h"
#include "gperftools/tcmalloc.h"

int main(int argc, char *argv[]) {
  DBAdaptor *dba;
  DNAAlignFeatureAdaptor *dafa;
  Slice *slice;
  Vector *features;
  int i;
  int failed;
  int failedTests = 0;
  
  initEnsC(argc, argv);

  dba = Test_initROEnsDB();

  slice = Test_getStandardSlice(dba);

  failedTests += ok(1, slice!=NULL);

  dafa = DBAdaptor_getDNAAlignFeatureAdaptor(dba);
  SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(dba);

  failedTests += ok(2, dafa!=NULL);

  //features =  Slice_getAllDNAAlignFeatures(slice,NULL,NULL, NULL,NULL);

  //Slice *slice3 = SliceAdaptor_fetchByRegion(sa,"chromosome","1",2,260000000,1,NULL,0);
  //Slice *slice2 = SliceAdaptor_fetchByRegion(sa,"chromosome","Y",1000000,4000000,1,NULL,0);
  features =  Slice_getAllDNAAlignFeatures(slice,NULL,NULL, NULL,NULL);

  failedTests += ok(3, features!=NULL);
  failedTests += ok(4, Vector_getNumElement(features)!=0);
//  MallocExtension_ReleaseFreeMemory();
//  tc_malloc_stats();
//  exit(1);


  failed = 0;

  for (i=0;i<Vector_getNumElement(features) && !failed;i++) {
    DNAAlignFeature *daf = Vector_getElementAt(features,i);
    Vector *ungapped;
    fprintf(stderr, "slice start %ld end %ld analysis %s\n", DNAAlignFeature_getStart(daf), DNAAlignFeature_getEnd(daf), Analysis_getLogicName(DNAAlignFeature_getAnalysis(daf)));
/*
    char *oldCigar = DNAAlignFeature_getCigarString(daf);

    //printf(" cigar = %s pre ungapped\n",DNAAlignFeature_getCigarString(daf));
    ungapped = DNAAlignFeature_getUngappedFeatures((BaseAlignFeature *)daf);
    if (!ungapped) failed = 1;
    //printf(" cigar = %s num ungapped %d\n",DNAAlignFeature_getCigarString(daf), Vector_getNumElement(ungapped));

    BaseAlignFeature_parseFeatures(daf,ungapped); 
    //printf(" cigar now = %s\n",DNAAlignFeature_getCigarString(daf));
// NIY Check that free func has been set
    Vector_setFreeFunc(ungapped, Object_freeImpl);
    Vector_free(ungapped);
    if (strcmp(oldCigar,DNAAlignFeature_getCigarString(daf))) {
      printf(" cigars different %s %s\n",oldCigar, DNAAlignFeature_getCigarString(daf));
      failed = 1;
    }
*/
  }
  failedTests += ok(5, !failed);

exit(1);

//  fprintf(stderr,"Clearing dafa cache\n");
//  DNAAlignFeatureAdaptor_clearCache(dafa);

//  Slice *slice2 = SliceAdaptor_fetchByRegion(sa,"chromosome","Y",1000000,4000000,1,NULL,0);
//  features =  Slice_getAllDNAAlignFeatures(slice2,NULL,NULL, NULL,NULL);
//  Slice *slice3 = SliceAdaptor_fetchByRegion(sa,"chromosome","1",10000000,100000000,1,NULL,0);
//  features =  Slice_getAllDNAAlignFeatures(slice3,NULL,NULL, NULL,NULL);
//  features =  Slice_getAllDNAAlignFeatures(slice,NULL,NULL, NULL,NULL);

  fprintf(stderr," \n\n\nTest 6 needs work - using SeqFeature transform and transfer methods - need to implement BaseAlignFeature ones\n\n\n");
  failed = 0;
  for (i=0;i<Vector_getNumElement(features) && !failed;i++) {
    DNAAlignFeature *daf = Vector_getElementAt(features,i);
    int start = DNAAlignFeature_getStart(daf);
    int end   = DNAAlignFeature_getEnd(daf);
    DNAAlignFeature *rdaf;

    // printf("slice start = %d end = %d id " IDFMTSTR "\n",start,end,DNAAlignFeature_getDbID(daf));
    // Temporary:
    rdaf = SeqFeature_transform(daf,"contig", NULL, NULL);
    // SHOULD BE THISrdaf = DNAAlignFeature_transform(daf,"contig", NULL, NULL);
    if (rdaf == NULL) {
//      printf("Feature didn't map\n");
    } else {

      // printf("rc id " IDFMTSTR " rc start = %d end = %d\n",BaseContig_getDbID(DNAAlignFeature_getContig(rdaf)),
      //       DNAAlignFeature_getStart(rdaf),DNAAlignFeature_getEnd(rdaf));
  
      // Temporary:
      daf = SeqFeature_transfer(rdaf, slice);
      // SHOULD BE daf = DNAAlignFeature_transfer(rdaf, slice);
      if (DNAAlignFeature_getStart(daf) != start ||
          DNAAlignFeature_getEnd(daf) != end) {
        // printf("slice start now = %d end = %d\n",DNAAlignFeature_getStart(daf),DNAAlignFeature_getEnd(daf));
        printf("Remapping to slice produced different coords for " IDFMTSTR "\n", DNAAlignFeature_getDbID(daf));
        failed =1;
      }
    }
  }
  failedTests += ok(6, !failed);
  return failedTests;
}
