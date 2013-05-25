#include <stdio.h>

#include "SliceAdaptor.h"
#include "ProteinAlignFeatureAdaptor.h"
#include "DBAdaptor.h"
#include "EnsC.h"
#include "DNAPepAlignFeature.h"

#include "BaseRODBTest.h"
#include "BaseRWDBTest.h"

int main(int argc, char *argv[]) {
  DBAdaptor *dba;
  DBAdaptor *writeDba;
  ProteinAlignFeatureAdaptor *pafa;
  Slice *slice;
  Vector *features;
  int i;
  int failed;
  
  initEnsC(argc, argv);

  dba = Test_initROEnsDB();

  writeDba = Test_initRWEnsDB();

  slice = Test_getStandardSlice(dba);

  ok(1, slice!=NULL);

  pafa = DBAdaptor_getProteinAlignFeatureAdaptor(writeDba);
  SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(dba);

  ok(2, pafa!=NULL);

  //features =  Slice_getAllDNAPepAlignFeatures(slice,NULL,NULL, NULL,NULL);

  //Slice *slice3 = SliceAdaptor_fetchByRegion(sa,"chromosome","1",2,260000000,1,NULL,0);
  Slice *slice2 = SliceAdaptor_fetchByRegion(sa,"chromosome","1",1000000,4000000,1,NULL,0);
  features =  Slice_getAllProteinAlignFeatures(slice2,NULL,NULL, NULL,NULL);

  ok(3, features!=NULL);
  ok(4, Vector_getNumElement(features)!=0);

  ProteinAlignFeatureAdaptor_store(pafa, features);

  return 0;
}
