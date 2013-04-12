#include <stdio.h>

#include "SliceAdaptor.h"
#include "DBAdaptor.h"
#include "EnsC.h"
#include "CoordSystemAdaptor.h"
#include "CoordSystem.h"

#include "BaseRODBTest.h"

int main(int argc, char *argv[]) {
  DBAdaptor *dba;
  CoordSystemAdaptor *csa;
  Vector *path;
  CoordSystem *cs;
  int i;
  int failed;
  
  initEnsC();

  dba = Test_initROEnsDB();

  csa = DBAdaptor_getCoordSystemAdaptor(dba);

  ok(1, csa!=NULL);

  Vector *csVec = CoordSystemAdaptor_fetchAllByName(csa,"chromosome");

  ok(2, Vector_getNumElement(csVec) > 0);

  cs = NULL;
  cs = CoordSystemAdaptor_fetchByName(csa,"chromosome",NULL);
  ok(3, cs != NULL);

  cs = CoordSystemAdaptor_fetchByName(csa,"notacoordsystemname",NULL);
  ok(4, cs == NULL);

  // Contigs should not have a version
  cs = CoordSystemAdaptor_fetchByName(csa,"contig","GRCh37");
  ok(5, cs == NULL);

  cs = CoordSystemAdaptor_fetchByName(csa,"contig",NULL);
  ok(6, cs != NULL);

  CoordSystem *cs2 = CoordSystemAdaptor_fetchByName(csa,"chromosome",NULL);
  path = CoordSystemAdaptor_getMappingPath(csa, cs,cs2);

  ok(7, Vector_getNumElement(path) > 0);

  return 0;
}
