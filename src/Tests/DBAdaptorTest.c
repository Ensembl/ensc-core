#include <stdio.h>

#include "DBAdaptor.h"
#include "EnsC.h"

#include "BaseTest.h"

int main(int argc, char *argv[]) {
  DBAdaptor *dba;

  initEnsC();

  dba = DBAdaptor_new("kaka.sanger.ac.uk","anonymous",NULL,"homo_sapiens_core_12_31",3306,NULL);

  ok(1,!strcmp("NCBI31",DBAdaptor_getAssemblyType(dba)));

  return 0;
}
