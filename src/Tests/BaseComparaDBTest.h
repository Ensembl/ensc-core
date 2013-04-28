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

Slice *Test_getStandardSlice(ComparaDBAdaptor *cdba) {
  Slice *slice;
  SliceAdaptor *sa;
  DBAdaptor *dba = ComparaDBAdaptor_getDBAdaptor(cdba,"Homo sapiens", "NCBI31");

  sa = DBAdaptor_getSliceAdaptor(dba);

  slice = SliceAdaptor_fetchByRegion(sa,"chromosome","1",500000,1000000,1,NULL,0);

  return slice;
}
#endif 
