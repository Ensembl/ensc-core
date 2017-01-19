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
#include "TopLevelAssemblyMapper.h"

#include "BaseRODBTest.h"


void printCoords(MapperRangeSet *results);
#define NumOutput(a) sizeof(a)/(sizeof(int)*4)

int main(int argc, char *argv[]) {
  DBAdaptor *dba;
  AssemblyMapperAdaptor *asma;
  int testNum = 1;

  initEnsC(argc, argv);

  dba = Test_initROEnsDB();

  //
  // 1 Test AssemblyMapperAdaptor constructor
  //
  asma = DBAdaptor_getAssemblyMapperAdaptor(dba);

  ok(testNum++, asma!=NULL);

  //
  // Test fetch_by_CoordSystems
  //

  CoordSystemAdaptor *csa = DBAdaptor_getCoordSystemAdaptor(dba);

  CoordSystemAdaptor_dumpCachedMappings(csa);

  CoordSystem *toplevelCs  = CoordSystemAdaptor_fetchByName(csa, "toplevel", NULL);
  CoordSystem *clnCs  = CoordSystemAdaptor_fetchByName(csa, "clone", NULL);
  CoordSystem *superctgCs = CoordSystemAdaptor_fetchByName(csa, "supercontig", NULL);

  TopLevelAssemblyMapper *clnToplevelMapper = (TopLevelAssemblyMapper *)AssemblyMapperAdaptor_fetchByCoordSystems(asma, toplevelCs, clnCs);
  TopLevelAssemblyMapper *superctgToplevelMapper = (TopLevelAssemblyMapper *)AssemblyMapperAdaptor_fetchByCoordSystems(asma, toplevelCs, superctgCs);

  ok(testNum++, clnToplevelMapper!=NULL); //  && $cln_toplevel_mapper->isa('Bio::EnsEMBL::TopLevelAssemblyMapper'));
  ok(testNum++, superctgToplevelMapper!=NULL); //  && $cln_toplevel_mapper->isa('Bio::EnsEMBL::TopLevelAssemblyMapper'));


//
// test db has chr 20  (50KB -> 62MB)
//

//
// Test map
//
  MapperRangeSet *coords = NULL;

  if (clnToplevelMapper) {
    fprintf(stderr, "MAP 'AL359765.6'->toplevel\n");
    coords = TopLevelAssemblyMapper_map(clnToplevelMapper,"AL359765.6", 1, 13780, 1, clnCs,  0, NULL);
    printCoords(coords);
    ok(testNum++, coords!=NULL);
  }

  if (superctgToplevelMapper) {
      fprintf(stderr, "MAP NT_028392->toplevel\n");
      coords = TopLevelAssemblyMapper_map(superctgToplevelMapper, "NT_028392", 600000, 1000000, 1, superctgCs, 0, NULL);
      printCoords(coords);
      ok(testNum++, coords!=NULL);
    }


//
// Test list_seq_regions
//
  Vector *seqRegions;
  int i;
    
  if (clnToplevelMapper) {
    seqRegions = TopLevelAssemblyMapper_listSeqRegions(clnToplevelMapper, "AL359765.6", 1, 13780, clnCs);
    ok(testNum++, seqRegions!=NULL && Vector_getNumElement(seqRegions) == 1 && !strcmp("20", Vector_getElementAt(seqRegions,0)));
    for (i=0;i<Vector_getNumElement(seqRegions); i++) {
      char *regionName = Vector_getElementAt(seqRegions, i);
      fprintf(stderr, "%s\n",regionName);
    }
  }

  if (superctgToplevelMapper) {
    seqRegions = TopLevelAssemblyMapper_listSeqRegions(superctgToplevelMapper, "NT_028392", 600000, 1000000, superctgCs);
    ok(testNum++, seqRegions!=NULL && Vector_getNumElement(seqRegions) == 1 && !strcmp("20", Vector_getElementAt(seqRegions,0)));
    for (i=0;i<Vector_getNumElement(seqRegions); i++) {
      char *regionName = Vector_getElementAt(seqRegions, i);
      fprintf(stderr, "%s\n",regionName);
    }
  }

//
// Test list_seq_ids
//
  Vector *ids;

  if (clnToplevelMapper) {
    ids = TopLevelAssemblyMapper_listIds(clnToplevelMapper, "AL359765.6", 1, 13780, clnCs);
    ok(testNum++, ids!=NULL && Vector_getNumElement(ids) == 1 && *((IDType *)Vector_getElementAt(ids,0)) == 469283 );
    for (i=0;i<Vector_getNumElement(ids); i++) {
      IDType id = *((IDType *)Vector_getElementAt(ids, i));
      fprintf(stderr, IDFMTSTR"\n",id);
    }
  }

  if (superctgToplevelMapper) {
    ids = TopLevelAssemblyMapper_listIds(superctgToplevelMapper, "NT_028392", 600000, 1000000, superctgCs);
    ok(testNum++, ids!=NULL && Vector_getNumElement(ids) == 1 && *((IDType *)Vector_getElementAt(ids,0)) == 469283 );
    for (i=0;i<Vector_getNumElement(ids); i++) {
      IDType id = *((IDType *)Vector_getElementAt(ids, i));
      fprintf(stderr, IDFMTSTR"\n",id);
    }
  }

// Test for a not implemented method
//  seqRegions = TopLevelAssemblyMapper_listContigIds(clnToplevelMapper, "AL359765.6", 1, 13780, 1);

  return 0;
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
