#include <stdio.h>

#include "SliceAdaptor.h"
#include "RepeatFeatureAdaptor.h"
#include "DBAdaptor.h"
#include "EnsC.h"
#include "RepeatFeature.h"

#include "BaseRODBTest.h"
#include "BaseRWDBTest.h"

int main(int argc, char *argv[]) {
  DBAdaptor *dba;
  DBAdaptor *writeDba;
  RepeatFeatureAdaptor *sfa;
  Slice *slice;
  Vector *features;
  int i;
  int failed;
  
  initEnsC(argc, argv);

  dba = Test_initROEnsDB();

  writeDba = Test_initRWEnsDB();

  slice = Test_getStandardSlice(dba);

  ok(1, slice!=NULL);

  sfa = DBAdaptor_getRepeatFeatureAdaptor(writeDba);
  SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(dba);

  ok(2, sfa!=NULL);

  //features =  Slice_getAllRepeatFeatures(slice,NULL,NULL, NULL,NULL);

  //Slice *slice3 = SliceAdaptor_fetchByRegion(sa,"chromosome","1",2,260000000,1,NULL,0);
  Slice *slice2 = SliceAdaptor_fetchByRegion(sa,"chromosome","1",1000000,4000000,1,NULL,0);
  features =  Slice_getAllRepeatFeatures(slice2,NULL,NULL, NULL);

  ok(3, features!=NULL);
  ok(4, Vector_getNumElement(features)!=0);

  RepeatFeatureAdaptor_store(sfa, features);

  return 0;
}
