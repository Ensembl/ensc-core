#include <stdio.h>

#include "DBAdaptor.h"
#include "EnsC.h"

#include "BaseTest.h"

int main(int argc, char *argv[]) {
  DBAdaptor *dba;

  initEnsC();

  dba = DBAdaptor_new("ensembldb.ensembl.org","anonymous",NULL,"homo_sapiens_core_29_35b",3306,NULL);

  ok(1,!strcmp("NCBI35",DBAdaptor_getAssemblyType(dba)));

  return 0;
}
