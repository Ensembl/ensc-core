#include <stdio.h>

#include "SliceAdaptor.h"
#include "DBAdaptor.h"
#include "EnsC.h"

#include "BaseTest.h"

int main(int argc, char *argv[]) {
  DBAdaptor *dba;
  Slice *slice;
  SliceAdaptor *sa;

  initEnsC();

  dba = DBAdaptor_new("kaka.sanger.ac.uk","anonymous",NULL,"homo_sapiens_core_12_31",3306,NULL);

  ok(1,!strcmp("NCBI31",DBAdaptor_getAssemblyType(dba)));

  sa = DBAdaptor_getSliceAdaptor(dba);

  ok(2, sa!=NULL);

  slice = SliceAdaptor_fetchByChrStartEnd(sa,"1",1,10000);

  ok(3, slice!=NULL);

  printf("slice name = %s\n",Slice_getName(slice));
  ok(4,!strcmp(Slice_getName(slice),"1.1-10000"));

  printf("slice seq = %s\n",Slice_getSeq(slice));


  return 0;
}
