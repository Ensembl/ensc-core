#ifndef __BASERODBTEST_H__
#define __BASERODBTEST_H__

#include "BaseTest.h"

#include "DBAdaptor.h"

DBAdaptor *Test_initROEnsDB() {
  DBAdaptor *dba;

//  dba = DBAdaptor_new("kaka.sanger.ac.uk","anonymous",NULL,"homo_sapiens_core_14_31",3306,NULL);
//  dba = DBAdaptor_new("ecs2b.internal.sanger.ac.uk","ensro",NULL,"human_NCBI33_raw",3306,NULL);
//  dba = DBAdaptor_new("127.0.0.1","ensro",NULL,"human_NCBI33_raw",13302,NULL);
//  dba = DBAdaptor_new("127.0.0.1","root",NULL,"test_core",3306,NULL);
//  dba = DBAdaptor_new("127.0.0.1","root",NULL,"homo_sapiens_core_18_34",3306,NULL);
//  dba = DBAdaptor_new("127.0.0.1","root",NULL,"steve_feb04_comp",3306,NULL);
  dba = DBAdaptor_new("ens-staging.internal.sanger.ac.uk","ensadmin","ensembl","homo_sapiens_core_71_37",3306,NULL);

  return dba;
}

Slice *Test_getStandardSlice(DBAdaptor *dba) {
  Slice *slice;
  SliceAdaptor *sa;

  sa = DBAdaptor_getSliceAdaptor(dba);

  slice = SliceAdaptor_fetchByChrStartEnd(sa,"2",1,100000000);

  return slice;
}
#endif 
