#include <stdio.h>

#include "SliceAdaptor.h"
#include "DBAdaptor.h"
#include "EnsC.h"
#include "SimpleFeature.h"

#include "BaseRODBTest.h"

int main(int argc, char *argv[]) {
  DBAdaptor *dba;
  SimpleFeatureAdaptor *sfa;
  Slice *slice;
  Vector *features;
  int i;
  int failed;
  
  initEnsC();

  dba = Test_initROEnsDB();

  slice = Test_getStandardSlice(dba);

  ok(1, slice!=NULL);

  sfa = DBAdaptor_getSimpleFeatureAdaptor(dba);

  ok(2, sfa!=NULL);

  features =  Slice_getAllSimpleFeatures(slice,"",NULL);

  ok(3, features!=NULL);
  ok(4, Vector_getNumElement(features)!=0);

  failed = 0;
  for (i=0;i<Vector_getNumElement(features) && !failed;i++) {
    SimpleFeature *sf = Vector_getElementAt(features,i);
    int start = SimpleFeature_getStart(sf);
    int end   = SimpleFeature_getEnd(sf);
    Vector *rsfVector;
    SimpleFeature *rsf;

    //printf("slice start = %d end = %d\n",start,end);
    rsfVector = SimpleFeature_transformToRawContig(sf);
    if (Vector_getNumElement(rsfVector) > 1) {
      printf("Feature mapped to more than one rawcontig\n");
      failed=1;
    }
    rsf = Vector_getElementAt(rsfVector,0);

    //printf("rc start = %d end = %d\n",SimpleFeature_getStart(rsf),SimpleFeature_getEnd(rsf));
    sf = SimpleFeature_transformToSlice(rsf, slice);
    if (SimpleFeature_getStart(sf) != start ||
        SimpleFeature_getEnd(sf) != end) {
      printf("Remapping to slice produced different coords\n");
      failed =1;
    }
  }
  ok(5, !failed);
  return 0;
}
