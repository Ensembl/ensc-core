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
  int testNum = 0;
  int testResult = 0;
  
  initEnsC(argc, argv);

  dba = Test_initROEnsDB();

  //
  // 1 Test AssemblyMapperAdaptor constructor
  //
  asma = DBAdaptor_getAssemblyMapperAdaptor(dba);

  testResult += ok(++testNum, asma!=NULL);


  //
  // 2 Test fetch_by_CoordSystems
  //
  CoordSystemAdaptor *csa = DBAdaptor_getCoordSystemAdaptor(dba);

  CoordSystem *ctgCs  = CoordSystemAdaptor_fetchByName(csa, "contig", NULL);
  CoordSystem *chrCs  = CoordSystemAdaptor_fetchByName(csa, "chromosome", NULL);
  CoordSystem *clnCs  = CoordSystemAdaptor_fetchByName(csa, "clone", NULL);
  MapperRangeSet *coords = NULL ;


/*
  testResult += ok(++testNum, $coords[0]->id() == 469271 );  #seq_region_id not name now.
  testResult += ok(++testNum, $coords[0]->start() == 91 );
  testResult += ok(++testNum, $coords[0]->end() == 200 );
  testResult += ok(++testNum, $coords[0]->strand() == 1 );

  testResult += ok(++testNum, $coords[1]->isa( "Bio::EnsEMBL::Mapper::Gap" ) );

  testResult += ok(++testNum, $coords[2]->id() == 469271);
  testResult += ok(++testNum, $coords[2]->start() == 201 );
  testResult += ok(++testNum, $coords[2]->end() == 400 );
  testResult += ok(++testNum, $coords[2]->strand() == -1 );

  testResult += ok(++testNum, $coords[4]->id() == 469282);
  testResult += ok(++testNum, $coords[4]->start() == 1 );
  testResult += ok(++testNum, $coords[4]->end() == 100 );
  testResult += ok(++testNum, $coords[4]->strand() == -1 );
*/
  
  
  AssemblyMapper *asmMapper = AssemblyMapperAdaptor_fetchByCoordSystems(asma, ctgCs, chrCs);

  testResult += ok(++testNum, asmMapper!=NULL); // && $asmMapper->isa('Bio::EnsEMBL::AssemblyMapper'));

  if (asmMapper)
    {
  
      // 
      // test db has chr 20  (50KB -> 62MB)
      //
  
      //
      // Test map
      //
  
      coords = AssemblyMapper_map(asmMapper, "20", 500001, 60000000, 1, chrCs, 0, NULL);
      testResult += ok(++testNum, coords!=NULL);
      printCoords(coords);
  
  
      coords = AssemblyMapper_map(asmMapper, "AL359765.6.1.13780", 1, 13780, 1, ctgCs, 0, NULL);
      testResult += ok(++testNum, coords!=NULL);
      printCoords(coords);
  
      AssemblyMapper *clnMapper = AssemblyMapperAdaptor_fetchByCoordSystems(asma, clnCs, ctgCs);
      coords = AssemblyMapper_map(clnMapper, "AL359765.6", 1, 20000, 1, clnCs, 0, NULL);
      testResult += ok(++testNum, coords!=NULL);
      printCoords(coords);
  
  
      // 
      // Test list_seq_regions
      //
  
      Vector *seqRegions = AssemblyMapper_listSeqRegions(asmMapper, "20", 500001, 60000000, chrCs);

      testResult += ok(++testNum, seqRegions != NULL);
      int i;
      for (i=0;i<Vector_getNumElement(seqRegions); i++) {
        char *regionName = Vector_getElementAt(seqRegions, i);
        fprintf(stderr, "%s\n",regionName);
      }
  
      seqRegions = AssemblyMapper_listSeqRegions(asmMapper, "AL359765.6", 1, 13780, ctgCs);
      testResult += ok(++testNum, seqRegions!=NULL);

      if (seqRegions) {
        for (i=0;i<Vector_getNumElement(seqRegions); i++) {
          char *regionName = Vector_getElementAt(seqRegions, i);
          fprintf(stderr, "%s\n",regionName);
        }
      }

  
      //
      // Test list_seq_ids
      //
  
      Vector *seqIds = AssemblyMapper_listIds(asmMapper, "20", 500001, 60000000, chrCs);

      testResult += ok(++testNum, seqIds!=NULL);

      if (seqIds) {
        for (i=0;i<Vector_getNumElement(seqIds); i++) {
          IDType regionId = *((IDType *)Vector_getElementAt(seqIds, i));
          fprintf(stderr, IDFMTSTR"\n",regionId);
        }
      }
  
      seqIds = AssemblyMapper_listIds(asmMapper, "AL359765.6", 1, 13780, ctgCs);
      testResult += ok(++testNum, seqIds!=NULL);

      if (seqIds) {
        for (i=0;i<Vector_getNumElement(seqIds); i++) {
          IDType regionId = *((IDType *)Vector_getElementAt(seqIds, i));
          fprintf(stderr, IDFMTSTR"\n",regionId);
        }
      }
    }
  
  return testResult;
}


void printCoords(MapperRangeSet *results) {
  int i;
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

int compareTransform(MapperRangeSet *results, int dest[][4], int nDest ) {
  int diff = 0;
  fprintf(stderr, "Number of results = %d nDest = %d\n",MapperRangeSet_getNumRange(results), nDest);
  if (MapperRangeSet_getNumRange(results) != nDest) {
    diff =1;

  } else{
    int i;
    for (i=0;i<MapperRangeSet_getNumRange(results) && !diff;i++) {
      MapperRange *range = MapperRangeSet_getRangeAt(results, i);

      switch (range->rangeType) {
        case MAPPERRANGE_COORD :
          {
            MapperCoordinate *mc = (MapperCoordinate *)range;
            fprintf(stderr, "Coord: "IDFMTSTR" %ld %ld %d\n", mc->id, mc->start, mc->end, mc->strand);
            if (dest[i][0] != mc->id || dest[i][1] != mc->start || dest[i][2] != mc->end || dest[i][3] != mc->strand) {
              diff=1;
            }
          }
          break;
        case MAPPERRANGE_GAP :
          {
            MapperGap *mg = (MapperGap *)range;

            fprintf(stderr, "Gap: %ld %ld\n", mg->start, mg->end);
            if (dest[i][0] != 0 || dest[i][1] != mg->start || dest[i][2] != mg->end || dest[i][3] != 0) {
              diff=1;
            }
          }
          break;
        default:
          {
            fprintf(stderr, "Unhandled range type %d\n",range->rangeType);
            diff=1;
          }
          break;
      }
    }
  }
  if (diff) {
    fprintf(stderr, "DIFFERENCE\n");
  }
  fprintf(stderr, "\n");
  return !diff;
}

