#include <stdio.h>

#include "SliceAdaptor.h"
#include "DBAdaptor.h"
#include "EnsC.h"
#include "SimpleFeature.h"
#include "BaseFeatureAdaptor.h"

#include "BaseRODBTest.h"

int main(int argc, char *argv[]) {
  DBAdaptor *dba;
  SimpleFeatureAdaptor *sfa;
  Slice *slice;
  Vector *features;
  int i;
  int failed;
  
  initEnsC(argc, argv);

  dba = Test_initROEnsDB();

  slice = Test_getStandardSlice(dba);

  ok(1, slice!=NULL);

  sfa = DBAdaptor_getSimpleFeatureAdaptor(dba);

  ok(2, sfa!=NULL);

  features =  Slice_getAllSimpleFeatures(slice, NULL, NULL, NULL);

  ok(3, features!=NULL);
  ok(4, Vector_getNumElement(features)!=0);

  failed = 0;
  for (i=0;i<Vector_getNumElement(features) && !failed;i++) {
    SimpleFeature *sf = Vector_getElementAt(features,i);
    long start = SimpleFeature_getStart(sf);
    long end   = SimpleFeature_getEnd(sf);
    Vector *rsfVector;
    SimpleFeature *rsf;

    //printf("slice start = %d end = %d\n",start,end);
    rsf = SeqFeature_transform(sf,"contig",NULL,NULL);

    if (rsf) {
//      printf("rc start = %ld end = %ld\n",SimpleFeature_getStart(rsf),SimpleFeature_getEnd(rsf));
    } else {
//      printf("no mapped feature\n");
    }

    if (rsf) {
      //sf = SeqFeature_transform(rsf,"chromosome",NULL,slice);
      sf = SeqFeature_transfer(rsf, slice);
      if (SimpleFeature_getStart(sf) != start ||
          SimpleFeature_getEnd(sf) != end) {
        printf("Remapping to slice produced different coords start %ld v %ld   end %ld v %ld\n", 
               SimpleFeature_getStart(sf),start, SimpleFeature_getEnd(sf),end);
        failed =1;
      }
    }
  }
  ok(5, !failed);
  return 0;
}
