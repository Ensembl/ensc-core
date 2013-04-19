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
  
  initEnsC();

  dba = Test_initROEnsDB();

  //
  // 1 Test AssemblyMapperAdaptor constructor
  //
  asma = DBAdaptor_getAssemblyMapperAdaptor(dba);

  ok(testNum++, asma!=NULL);


  //
  // 2 Test fetch_by_CoordSystems
  //

  CoordSystemAdaptor *csa = DBAdaptor_getCoordSystemAdaptor(dba);

  CoordSystemAdaptor_dumpCachedMappings(csa);

  CoordSystem *chrCs  = CoordSystemAdaptor_fetchByName(csa, "chromosome", NULL);
  CoordSystem *clnCs  = CoordSystemAdaptor_fetchByName(csa, "clone", NULL);
  CoordSystem *sCtgCs = CoordSystemAdaptor_fetchByName(csa, "supercontig", NULL);

  ChainedAssemblyMapper *asmMapper =  (ChainedAssemblyMapper *)AssemblyMapperAdaptor_fetchByCoordSystems(asma, clnCs, chrCs);

  ok(testNum++,  asmMapper!=NULL); // Need to make it an object before can do this && asmMapper->objectType == ( "Bio::EnsEMBL::ChainedAssemblyMapper" ));
  
  ChainedAssemblyMapper *chrSCtgMapper = (ChainedAssemblyMapper *)AssemblyMapperAdaptor_fetchByCoordSystems(asma, chrCs, sCtgCs);

  ok(testNum++, chrSCtgMapper!=NULL);// && $chr_sctg_mapper->isa('Bio::EnsEMBL::ChainedAssemblyMapper'));

//
// test db has chr 20  (50KB -> 62MB)
//

  MapperRangeSet *coords;

  fprintf(stderr,"MAP 20->clone\n");
  coords = ChainedAssemblyMapper_map(asmMapper, "20", 500001, 60000000, 1, chrCs, 0, NULL);
  ok(testNum++, coords!=NULL);
  printCoords(coords);


  fprintf(stderr,"MAP 'AL359765.6'->chromosome\n");
  coords = ChainedAssemblyMapper_map(asmMapper, "AL359765.6", 1, 13780, 1, clnCs, 0, NULL);
  ok(testNum++, coords!=NULL);
  printCoords(coords);

  fprintf(stderr,"MAP 20->supercontig\n");
  coords = ChainedAssemblyMapper_map(chrSCtgMapper, "20", 500001, 60000000, 1, chrCs, 0, NULL);
  ok(testNum++, coords!=NULL);
  printCoords(coords);


  //
  // Test list_seq_regions
  //
  fprintf(stderr,"Starting list tests\n");

  Vector *seqRegions = ChainedAssemblyMapper_listSeqRegions(asmMapper, "20", 500001, 60000000, chrCs);
  ok(testNum++, seqRegions != NULL);
  int i;
  for (i=0;i<Vector_getNumElement(seqRegions); i++) {
    char *regionName = Vector_getElementAt(seqRegions, i);
    fprintf(stderr, "%s\n",regionName);
  }

  seqRegions = ChainedAssemblyMapper_listSeqRegions(asmMapper, "AL359765.6", 1, 13780, clnCs);
  ok(testNum++, seqRegions!=NULL);
  for (i=0;i<Vector_getNumElement(seqRegions); i++) {
    char *regionName = Vector_getElementAt(seqRegions, i);
    fprintf(stderr, "%s\n",regionName);
  }


  seqRegions = ChainedAssemblyMapper_listSeqRegions(chrSCtgMapper, "NT_028392", 600000, 1000000, sCtgCs);
  ok(testNum++, seqRegions!=NULL);
  for (i=0;i<Vector_getNumElement(seqRegions); i++) {
    char *regionName = Vector_getElementAt(seqRegions, i);
    fprintf(stderr, "%s\n",regionName);
  }

  seqRegions = ChainedAssemblyMapper_listSeqRegions(chrSCtgMapper, "20", 3000000, 31000000, chrCs);
  ok(testNum++, seqRegions!=NULL);
  for (i=0;i<Vector_getNumElement(seqRegions); i++) {
    char *regionName = Vector_getElementAt(seqRegions, i);
    fprintf(stderr, "%s\n",regionName);
  }



  //
  // Test list_seq_ids
  //

  Vector *seqIds = ChainedAssemblyMapper_listIds(asmMapper, "20", 500001, 60000000, chrCs);

  ok(testNum++, seqIds!=NULL);
  for (i=0;i<Vector_getNumElement(seqIds); i++) {
    IDType regionId = *((IDType *)Vector_getElementAt(seqIds, i));
    fprintf(stderr, IDFMTSTR"\n",regionId);
  }

  seqIds = ChainedAssemblyMapper_listIds(asmMapper, "AL359765.6", 1, 13780, clnCs);
  ok(testNum++, seqIds!=NULL);
  for (i=0;i<Vector_getNumElement(seqIds); i++) {
    IDType regionId = *((IDType *)Vector_getElementAt(seqIds, i));
    fprintf(stderr, IDFMTSTR"\n",regionId);
  }

  seqIds = ChainedAssemblyMapper_listIds(chrSCtgMapper, "NT_028392", 600000, 1000000, sCtgCs);
  ok(testNum++, seqIds!=NULL);
  for (i=0;i<Vector_getNumElement(seqIds); i++) {
    IDType regionId = *((IDType *)Vector_getElementAt(seqIds, i));
    fprintf(stderr, IDFMTSTR"\n",regionId);
  }

  seqIds = ChainedAssemblyMapper_listIds(chrSCtgMapper, "20", 3000000, 31000000, chrCs);
  ok(testNum++, seqIds!=NULL);
  for (i=0;i<Vector_getNumElement(seqIds); i++) {
    IDType regionId = *((IDType *)Vector_getElementAt(seqIds, i));
    fprintf(stderr, IDFMTSTR"\n",regionId);
  }

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


