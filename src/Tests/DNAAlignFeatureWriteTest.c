#include <stdio.h>

#include "SliceAdaptor.h"
#include "DNAAlignFeatureAdaptor.h"
#include "DBAdaptor.h"
#include "EnsC.h"
#include "DNAAlignFeature.h"

#include "BaseRODBTest.h"
#include "BaseRWDBTest.h"
#include "gperftools/tcmalloc.h"

int main(int argc, char *argv[]) {
  DBAdaptor *dba;
  DBAdaptor *writeDba;
  DNAAlignFeatureAdaptor *dafa;
  Slice *slice;
  Vector *features;
  int i;
  int failed;
  
  initEnsC(argc, argv);

  dba = Test_initROEnsDB();

  writeDba = Test_initRWEnsDB();

  slice = Test_getStandardSlice(dba);

  ok(1, slice!=NULL);

  dafa = DBAdaptor_getDNAAlignFeatureAdaptor(writeDba);
  SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(dba);

  ok(2, dafa!=NULL);

  //features =  Slice_getAllDNAAlignFeatures(slice,NULL,NULL, NULL,NULL);

  //Slice *slice3 = SliceAdaptor_fetchByRegion(sa,"chromosome","1",2,260000000,1,NULL,0);
  Slice *slice2 = SliceAdaptor_fetchByRegion(sa,"chromosome","1",1000000,4000000,1,NULL,0);
  features =  Slice_getAllDNAAlignFeatures(slice2,NULL,NULL, NULL,NULL);

  ok(3, features!=NULL);
  ok(4, Vector_getNumElement(features)!=0);

  DNAAlignFeatureAdaptor_store(dafa, features);

  return 0;
}
