#ifndef __BASERWDBTEST_H__
#define __BASERWDBTEST_H__

#include "BaseTest.h"

#include "DBAdaptor.h"

DBAdaptor *Test_initRWEnsDB() {
  DBAdaptor *dba;

//  dba = DBAdaptor_new("kaka.sanger.ac.uk","anonymous",NULL,"homo_sapiens_core_14_31",3306,NULL);
//  dba = DBAdaptor_new("ecs2b.internal.sanger.ac.uk","ensro",NULL,"human_NCBI33_raw",3306,NULL);
//  dba = DBAdaptor_new("127.0.0.1","ensro",NULL,"human_NCBI33_raw",13302,NULL);
  dba = DBAdaptor_new("127.0.0.1","root",NULL,"test_core",3306,NULL);

  return dba;
}

Slice *Test_getStandardSlice(DBAdaptor *dba) {
  Slice *slice;
  SliceAdaptor *sa;

  sa = DBAdaptor_getSliceAdaptor(dba);

  slice = SliceAdaptor_fetchByChrStartEnd(sa,"1",1,5000000);

  return slice;
}
#endif 
