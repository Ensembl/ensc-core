#include <stdio.h>

#include "SliceAdaptor.h"
#include "DBAdaptor.h"
#include "EnsC.h"
#include "RepeatFeature.h"

#include "BaseRODBTest.h"

int main(int argc, char *argv[]) {
  DBAdaptor *dba;
  RepeatFeatureAdaptor *rfa;
  Slice *slice;
  Vector *features;
  int i;
  int failed;
  
  initEnsC();

  dba = Test_initROEnsDB();

  slice = Test_getStandardSlice(dba);

  ok(1, slice!=NULL);

  rfa = DBAdaptor_getRepeatFeatureAdaptor(dba);

  ok(2, rfa!=NULL);

  features =  Slice_getAllRepeatFeatures(slice,"");

  ok(3, features!=NULL);
  ok(4, Vector_getNumElement(features)!=0);

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
  ok(5, !failed);
  return 0;
}
