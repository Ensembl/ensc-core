/*
 * Copyright [1999-2024] EMBL-European Bioinformatics Institute
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *      http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

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
  int testResult = 0;
  
  initEnsC(argc, argv);

  dba = Test_initROEnsDB();

  csa = DBAdaptor_getCoordSystemAdaptor(dba);

  testResult += ok(1, csa!=NULL);

  Vector *csVec = CoordSystemAdaptor_fetchAllByName(csa,"chromosome");

  testResult += ok(2, Vector_getNumElement(csVec) > 0);

  cs = NULL;
  cs = CoordSystemAdaptor_fetchByName(csa,"chromosome",NULL);
  testResult += ok(3, cs != NULL);

  cs = CoordSystemAdaptor_fetchByName(csa,"notacoordsystemname",NULL);
  testResult += ok(4, cs == NULL);

  // Contigs should not have a version
  cs = CoordSystemAdaptor_fetchByName(csa,"contig","GRCh37");
  testResult += ok(5, cs == NULL);

  cs = CoordSystemAdaptor_fetchByName(csa,"contig",NULL);
  testResult += ok(6, cs != NULL);

  CoordSystem *cs2 = CoordSystemAdaptor_fetchByName(csa,"chromosome",NULL);
  path = CoordSystemAdaptor_getMappingPath(csa, cs,cs2);

  testResult += ok(7, Vector_getNumElement(path) > 0);

  return testResult;
}
