#ifndef __BASECOMPARADBTEST_H__
#define __BASECOMPARADBTEST_H__

#include "BaseTest.h"

#include "ComparaDBAdaptor.h"
#include "DBAdaptor.h"

ComparaDBAdaptor *Test_initComparaDB() {
  ComparaDBAdaptor *cdba;

  cdba = ComparaDBAdaptor_new("kaka.sanger.ac.uk","anonymous",NULL,"ensembl_compara_14_1",3306,"Tests/Compara.conf");

  return cdba;
}

Slice *Test_getStandardSlice(DBAdaptor *dba) {
  Slice *slice;
  SliceAdaptor *sa;

  sa = DBAdaptor_getSliceAdaptor(dba);

  slice = SliceAdaptor_fetchByChrStartEnd(sa,"1",1,5000000);

  return slice;
}
#endif 
