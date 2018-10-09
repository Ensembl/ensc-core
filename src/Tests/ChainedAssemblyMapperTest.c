/*
 * Copyright [1999-2017] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
#include "AssemblyMapperAdaptor.h"
#include "CoordSystem.h"
#include "CoordSystemAdaptor.h"
#include "ChainedAssemblyMapper.h"

#include "BaseRODBTest.h"


int compareTransform(MapperRangeSet *results, int dest[][4], int nDest );
void printCoords(MapperRangeSet *results);
#define NumOutput(a) sizeof(a)/(sizeof(int)*4)

int main(int argc, char *argv[]) {
  DBAdaptor *dba;
  AssemblyMapperAdaptor *asma;
  int testNum = 1;
  int testResult = 0;
  
  initEnsC(argc, argv);

  dba = Test_initROEnsDB();

  //
  // 1 Test AssemblyMapperAdaptor constructor
  //
  asma = DBAdaptor_getAssemblyMapperAdaptor(dba);

  testResult += ok(testNum++, asma!=NULL);


  //
  // 2 Test fetch_by_CoordSystems
  //

  CoordSystemAdaptor *csa = DBAdaptor_getCoordSystemAdaptor(dba);

  CoordSystemAdaptor_dumpCachedMappings(csa);

  CoordSystem *chrCs  = CoordSystemAdaptor_fetchByName(csa, "chromosome", NULL);
  CoordSystem *clnCs  = CoordSystemAdaptor_fetchByName(csa, "clone", NULL);
  CoordSystem *sCtgCs = CoordSystemAdaptor_fetchByName(csa, "supercontig", NULL);

  ChainedAssemblyMapper *asmMapper =  (ChainedAssemblyMapper *)AssemblyMapperAdaptor_fetchByCoordSystems(asma, clnCs, chrCs);

  testResult += ok(testNum++,  asmMapper!=NULL); // Need to make it an object before can do this && asmMapper->objectType == ( "Bio::EnsEMBL::ChainedAssemblyMapper" ));
  
  ChainedAssemblyMapper *chrSCtgMapper = (ChainedAssemblyMapper *)AssemblyMapperAdaptor_fetchByCoordSystems(asma, chrCs, sCtgCs);

  testResult += ok(testNum++, chrSCtgMapper!=NULL);// && $chr_sctg_mapper->isa('Bio::EnsEMBL::ChainedAssemblyMapper'));

//
// test db has chr 20  (50KB -> 62MB)
//

  MapperRangeSet *coords;

  if (asmMapper)
    {
      fprintf(stderr,"MAP 20->clone\n");
      coords = ChainedAssemblyMapper_map(asmMapper, "20", 500001, 60000000, 1, chrCs, 0, NULL);
      testResult += ok(testNum++, coords!=NULL);
      printCoords(coords);
    }

  if (asmMapper)
    {
      fprintf(stderr,"MAP 'AL359765.6'->chromosome\n");
      coords = ChainedAssemblyMapper_map(asmMapper, "AL359765.6", 1, 13780, 1, clnCs, 0, NULL);
      testResult += ok(testNum++, coords!=NULL);
      printCoords(coords);
    }

  if (chrSCtgMapper)
    {
      fprintf(stderr,"MAP 20->supercontig\n");
      coords = ChainedAssemblyMapper_map(chrSCtgMapper, "20", 500001, 60000000, 1, chrCs, 0, NULL);
      testResult += ok(testNum++, coords!=NULL);
      printCoords(coords);
    }

  //
  // Test list_seq_regions
  //
  fprintf(stderr,"Starting list tests\n");
  int i;

  if (asmMapper)
    {
      Vector *seqRegions = ChainedAssemblyMapper_listSeqRegions(asmMapper, "20", 500001, 60000000, chrCs);
      testResult += ok(testNum++, seqRegions != NULL);
      for (i=0;i<Vector_getNumElement(seqRegions); i++) {
        char *regionName = Vector_getElementAt(seqRegions, i);
        fprintf(stderr, "%s\n",regionName);
      }
    }

  if (asmMapper)
    {
      Vector *seqRegions = ChainedAssemblyMapper_listSeqRegions(asmMapper, "AL359765.6", 1, 13780, clnCs);
      testResult += ok(testNum++, seqRegions!=NULL);
      for (i=0;i<Vector_getNumElement(seqRegions); i++) {
        char *regionName = Vector_getElementAt(seqRegions, i);
        fprintf(stderr, "%s\n",regionName);
      }
    }


  if (chrSCtgMapper)
    {
      Vector *seqRegions = ChainedAssemblyMapper_listSeqRegions(chrSCtgMapper, "NT_028392", 600000, 1000000, sCtgCs);
      testResult += ok(testNum++, seqRegions!=NULL);
      for (i=0;i<Vector_getNumElement(seqRegions); i++) {
        char *regionName = Vector_getElementAt(seqRegions, i);
        fprintf(stderr, "%s\n",regionName);
      }
    }

  if (chrSCtgMapper)
    {
      Vector *seqRegions = ChainedAssemblyMapper_listSeqRegions(chrSCtgMapper, "20", 3000000, 31000000, chrCs);
      testResult += ok(testNum++, seqRegions!=NULL);
      for (i=0;i<Vector_getNumElement(seqRegions); i++) {
        char *regionName = Vector_getElementAt(seqRegions, i);
        fprintf(stderr, "%s\n",regionName);
      }
    }



  //
  // Test list_seq_ids
  //

  if (asmMapper)
    {
      Vector *seqIds = ChainedAssemblyMapper_listIds(asmMapper, "20", 500001, 60000000, chrCs);

      testResult += ok(testNum++, seqIds!=NULL);
      for (i=0;i<Vector_getNumElement(seqIds); i++) {
        IDType regionId = *((IDType *)Vector_getElementAt(seqIds, i));
        fprintf(stderr, IDFMTSTR"\n",regionId);
      }
    }

  if (asmMapper)
    {
      Vector *seqIds = ChainedAssemblyMapper_listIds(asmMapper, "AL359765.6", 1, 13780, clnCs);
      testResult += ok(testNum++, seqIds!=NULL);
      for (i=0;i<Vector_getNumElement(seqIds); i++) {
        IDType regionId = *((IDType *)Vector_getElementAt(seqIds, i));
        fprintf(stderr, IDFMTSTR"\n",regionId);
      }
    }

  if (chrSCtgMapper)
    {
      Vector *seqIds = ChainedAssemblyMapper_listIds(chrSCtgMapper, "NT_028392", 600000, 1000000, sCtgCs);
      testResult += ok(testNum++, seqIds!=NULL);
      for (i=0;i<Vector_getNumElement(seqIds); i++) {
        IDType regionId = *((IDType *)Vector_getElementAt(seqIds, i));
        fprintf(stderr, IDFMTSTR"\n",regionId);
      }
    }

  if (chrSCtgMapper)
    {
      Vector *seqIds = ChainedAssemblyMapper_listIds(chrSCtgMapper, "20", 3000000, 31000000, chrCs);
      testResult += ok(testNum++, seqIds!=NULL);
      for (i=0;i<Vector_getNumElement(seqIds); i++) {
        IDType regionId = *((IDType *)Vector_getElementAt(seqIds, i));
        fprintf(stderr, IDFMTSTR"\n",regionId);
      }
    }

  return testResult;
  
}


void printCoords(MapperRangeSet *results) {
  int i;
  if (results) {
    for (i=0;i<MapperRangeSet_getNumRange(results);i++) {
      MapperRange *range = MapperRangeSet_getRangeAt(results, i);
      switch (range->rangeType) {
      case MAPPERRANGE_COORD :
        {
          MapperCoordinate *mc = (MapperCoordinate *)range;
          fprintf(stderr, "Coord: "IDFMTSTR" %ld %ld %d\n", mc->id, mc->start, mc->end, mc->strand);
        }
        break;
      case MAPPERRANGE_GAP :
        {
          MapperGap *mg = (MapperGap *)range;
          fprintf(stderr, "Gap: %ld %ld\n", mg->start, mg->end);
        }
        break;
      default:
        {
          fprintf(stderr, "Unhandled range type %d\n",range->rangeType);
        }
        break;
      }
    }
  }
}


