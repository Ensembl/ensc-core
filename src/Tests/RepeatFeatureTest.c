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
#include "RepeatFeature.h"

#include "BaseRODBTest.h"

int main(int argc, char *argv[]) {
  int testResult = 0;
  DBAdaptor *dba;
  RepeatFeatureAdaptor *rfa;
  Slice *slice;
  Vector *features;
  int i;
  int failed;
  
  initEnsC(argc, argv);

  dba = Test_initROEnsDB();

  slice = Test_getStandardSlice(dba);

  testResult += ok(1, slice!=NULL);

  rfa = DBAdaptor_getRepeatFeatureAdaptor(dba);

  testResult += ok(2, rfa!=NULL);

  features =  Slice_getAllRepeatFeatures(slice,NULL,NULL, NULL);

  testResult += ok(3, features!=NULL);
  testResult += ok(4, Vector_getNumElement(features)!=0);

  failed = 0;
  for (i=0;i<Vector_getNumElement(features) && !failed;i++) {
    RepeatFeature *rf = Vector_getElementAt(features,i);
    int start = RepeatFeature_getStart(rf);
    int end   = RepeatFeature_getEnd(rf);
    Vector *rrfVector;
    RepeatFeature *rrf;

    printf("slice start = %d end = %d\n",start,end);
/*
    rrfVector = RepeatFeature_transformToRawContig(rf);
    if (Vector_getNumElement(rrfVector) > 1) {
      printf("Feature mapped to more than one rawcontig\n");
      failed=1;
    }
    rrf = Vector_getElementAt(rrfVector,0);

    //printf("rc start = %d end = %d\n",RepeatFeature_getStart(rrf),RepeatFeature_getEnd(rrf));
    rf = RepeatFeature_transformToSlice(rrf, slice);
    if (RepeatFeature_getStart(rf) != start ||
        RepeatFeature_getEnd(rf) != end) {
      printf("Remapping to slice produced different coords\n");
      failed =1;
    }
*/
  }
  testResult += ok(5, !failed);
  return testResult;
}
