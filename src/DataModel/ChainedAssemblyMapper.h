#ifndef __CHAINEDASSEMBLYMAPPER_H__
#define __CHAINEDASSEMBLYMAPPER_H__

#include "DataModelTypes.h"
#include "EnsC.h"

#include "Vector.h"
#include "IDHash.h"
#include "Mapper.h"
#include "AdaptorTypes.h"
#include "MapperRangeSet.h"
#include "RangeRegistry.h"

#define MAXCHUNK 2000

struct ChainedAssemblyMapperStruct {
  AssemblyMapperAdaptor *adaptor;
  Mapper *firstMiddleMapper;
  Mapper *lastMiddleMapper;
  Mapper *firstLastMapper;
  RangeRegistry *firstRegister;
  RangeRegistry *lastRegister;
  CoordSystem *firstCoordSystem;
  CoordSystem *middleCoordSystem;
  CoordSystem *lastCoordSystem;
  int maxPairCount;
};

#define ChainedAssemblyMapper_setAdaptor(am, ad) (am)->adaptor = (ad)
#define ChainedAssemblyMapper_getAdaptor(am) (am)->adaptor

#define ChainedAssemblyMapper_setLastCoordSystem(am, cs) (am)->lastCoordSystem = (cs)
#define ChainedAssemblyMapper_getLastCoordSystem(am) (am)->lastCoordSystem

#define ChainedAssemblyMapper_setMiddleCoordSystem(am, cs) (am)->middleCoordSystem = (cs)
#define ChainedAssemblyMapper_getMiddleCoordSystem(am) (am)->middleCoordSystem

#define ChainedAssemblyMapper_setFirstCoordSystem(am, cs) (am)->firstCoordSystem = (cs)
#define ChainedAssemblyMapper_getFirstCoordSystem(am) (am)->firstCoordSystem

#define ChainedAssemblyMapper_setMaxPairCount(am, mp) (am)->maxPairCount = (mp)
#define ChainedAssemblyMapper_getMaxPairCount(am) (am)->maxPairCount

#define ChainedAssemblyMapper_setFirstLastMapper(am, m) (am)->firstLastMapper = (m)
#define ChainedAssemblyMapper_getFirstLastMapper(am) (am)->firstLastMapper

#define ChainedAssemblyMapper_setFirstMiddleMapper(am, m) (am)->firstMiddleMapper = (m)
#define ChainedAssemblyMapper_getFirstMiddleMapper(am) (am)->firstMiddleMapper

#define ChainedAssemblyMapper_setLastMiddleMapper(am, m) (am)->lastMiddleMapper = (m)
#define ChainedAssemblyMapper_getLastMiddleMapper(am) (am)->lastMiddleMapper

#define ChainedAssemblyMapper_setFirstRegistry(am, cr) (am)->firstRegister = (cr)
#define ChainedAssemblyMapper_getFirstRegistry(am) (am)->firstRegister

#define ChainedAssemblyMapper_setLastRegistry(am, cr) (am)->lastRegister = (cr)
#define ChainedAssemblyMapper_getLastRegistry(am) (am)->lastRegister

ChainedAssemblyMapper *ChainedAssemblyMapper_new(AssemblyMapperAdaptor *ama, Vector *coordSystems);
void ChainedAssemblyMapper_registerAll(ChainedAssemblyMapper *cam);
void ChainedAssemblyMapper_flush(ChainedAssemblyMapper *cam);
int ChainedAssemblyMapper_getSize(ChainedAssemblyMapper *cam);
MapperRangeSet *ChainedAssemblyMapper_map(ChainedAssemblyMapper *cam, char *frmSeqRegionName, long frmStart, long frmEnd, int frmStrand,
                              CoordSystem *frmCs, int fastmap, Slice *toSlice);
MapperRangeSet *ChainedAssemblyMapper_fastMap(ChainedAssemblyMapper *cam, char *frmSeqRegionName, long frmStart, long frmEnd, int frmStrand,
                              CoordSystem *frmCs);
Vector *ChainedAssemblyMapper_listIds(ChainedAssemblyMapper *cam, char *fromSeqRegionName, long frmStart, long frmEnd, CoordSystem *frmCs);
Vector *ChainedAssemblyMapper_listSeqRegions(ChainedAssemblyMapper *cam, char *fromSeqRegionName, long frmStart, long frmEnd, CoordSystem *frmCs);

#endif


