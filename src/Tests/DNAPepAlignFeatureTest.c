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
#include "DBAdaptor.h"
#include "EnsC.h"
#include "DNAPepAlignFeature.h"
#include "ProcUtil.h"

#include "BaseRODBTest.h"

int main(int argc, char *argv[]) {
  int failedTests = 0;
  DBAdaptor *dba;
  ProteinAlignFeatureAdaptor *pafa;
  Slice *slice;
  Vector *features;
  int i;
  int failed;
  
  initEnsC(argc, argv);

  dba = Test_initROEnsDB();

  slice = Test_getStandardSlice(dba);

  failedTests += ok(1, slice!=NULL);

  pafa = DBAdaptor_getProteinAlignFeatureAdaptor(dba);

  failedTests += ok(2, pafa!=NULL);

  features =  Slice_getAllProteinAlignFeatures(slice,NULL,NULL,NULL,NULL);

  failedTests += ok(3, features!=NULL);
  failedTests += ok(4, Vector_getNumElement(features)!=0);

  failed = 0;
  for (i=0;i<Vector_getNumElement(features) && !failed;i++) {
    DNAPepAlignFeature *daf = Vector_getElementAt(features,i);
    Vector *ungapped;
    char *oldCigar = DNAPepAlignFeature_getCigarString(daf);

    ungapped = DNAPepAlignFeature_getUngappedFeatures((BaseAlignFeature *)daf);
    if (!ungapped) failed = 1;
    // printf(" cigar = %s num ungapped %d\n",DNAPepAlignFeature_getCigarString(daf), Vector_getNumElement(ungapped));

    BaseAlignFeature_parseFeatures((BaseAlignFeature*)daf,ungapped); 
    // printf(" cigar now = %s\n",DNAPepAlignFeature_getCigarString(daf));
// NIY Make sure free func has been set for vector
    Vector_free(ungapped);
    if (strcmp(oldCigar,DNAPepAlignFeature_getCigarString(daf))) {
      printf(" cigars different %s %s\n",oldCigar, DNAPepAlignFeature_getCigarString(daf));
      failed = 1;
    }
  }
  failedTests += ok(5, !failed);

  failed = 0;
  for (i=0;i<Vector_getNumElement(features) && !failed;i++) {
    DNAPepAlignFeature *daf = Vector_getElementAt(features,i);
    int start = DNAPepAlignFeature_getStart(daf);
    int end   = DNAPepAlignFeature_getEnd(daf);
    Vector *rdafVector;
    DNAPepAlignFeature *rdaf;

    // printf("slice start = %d end = %d\n",start,end);
    // Temporary:
    rdaf = (DNAPepAlignFeature*)SeqFeature_transform((SeqFeature*)daf, "contig", NULL, NULL);
/*
    if (Vector_getNumElement(rdafVector) > 1) {
      printf("Feature mapped to more than one rawcontig\n");
      failed=1;
    }
    rdaf = Vector_getElementAt(rdafVector,0);
*/
    if (rdaf == NULL) {
      printf("Feature didn't map\n");
    } else {

    // printf("rc start = %d end = %d\n",DNAPepAlignFeature_getStart(rdaf),DNAPepAlignFeature_getEnd(rdaf));
      // Temporary:
      daf = (DNAPepAlignFeature*)SeqFeature_transfer((SeqFeature*)rdaf, slice);
      if (DNAPepAlignFeature_getStart(daf) != start ||
          DNAPepAlignFeature_getEnd(daf) != end) {
        printf("Remapping to slice produced different coords\n");
        failed =1;
      }
    }
  }
  failedTests += ok(6, !failed);

  Vector_free(features);

  ProcUtil_mallInfo();


  return failedTests;
}
