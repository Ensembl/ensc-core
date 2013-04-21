#ifndef __CHAINEDASSEMBLYMAPPER_H__
#define __CHAINEDASSEMBLYMAPPER_H__

#include "BaseAssemblyMapper.h"
#include "DataModelTypes.h"
#include "EnsC.h"

#include "Vector.h"
#include "IDHash.h"
#include "Mapper.h"
#include "AdaptorTypes.h"
#include "MapperRangeSet.h"
#include "RangeRegistry.h"

#define MAXCHUNK 2000

BASEASSEMBLYMAPPERFUNC_TYPES(ChainedAssemblyMapper)

typedef struct ChainedAssemblyMapperFuncsStruct {
  BASEASSEMBLYMAPPERFUNCS_DATA(ChainedAssemblyMapper)
} ChainedAssemblyMapperFuncs;


#define FUNCSTRUCTTYPE ChainedAssemblyMapperFuncs
struct ChainedAssemblyMapperStruct {
  BASEASSEMBLYMAPPER_DATA
  Mapper *firstLastMapper; // Deliberately put here so will be equivalent to mapper in AssemblyMapper
  Mapper *firstMiddleMapper;
  Mapper *lastMiddleMapper;
  RangeRegistry *firstRegister;
  RangeRegistry *lastRegister;
  CoordSystem *firstCoordSystem;
  CoordSystem *middleCoordSystem;
  CoordSystem *lastCoordSystem;
};
#undef FUNCSTRUCTTYPE

#define ChainedAssemblyMapper_getSeqRegionId(am, srName, cs) BaseAssemblyMapper_getSeqRegionId((BaseAssemblyMapper *)(am), (srName), (cs))

#define ChainedAssemblyMapper_setAdaptor(am, ad) BaseAssemblyMapper_setAdaptor(am,ad)
#define ChainedAssemblyMapper_getAdaptor(am) BaseAssemblyMapper_getAdaptor(am)

#define ChainedAssemblyMapper_setLastCoordSystem(am, cs) (am)->lastCoordSystem = (cs)
#define ChainedAssemblyMapper_getLastCoordSystem(am) (am)->lastCoordSystem

#define ChainedAssemblyMapper_setMiddleCoordSystem(am, cs) (am)->middleCoordSystem = (cs)
#define ChainedAssemblyMapper_getMiddleCoordSystem(am) (am)->middleCoordSystem

#define ChainedAssemblyMapper_setFirstCoordSystem(am, cs) (am)->firstCoordSystem = (cs)
#define ChainedAssemblyMapper_getFirstCoordSystem(am) (am)->firstCoordSystem

// perl assumes first cs = assembled cs
#define ChainedAssemblyMapper_getAssembledCoordSystem(am) (am)->firstCoordSystem

// perl assumes last cs = component cs
#define ChainedAssemblyMapper_getComponentCoordSystem(am) (am)->lastCoordSystem

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

#define ChainedAssemblyMapper_registerAll(am) ((am)->funcs->registerAll((am)))

#define ChainedAssemblyMapper_flush(am) ((am)->funcs->flush((am)))

#define ChainedAssemblyMapper_free(am) ((am)->funcs->free((am)))

#define ChainedAssemblyMapper_map(am, fsrName, fStart, fEnd, fStrand, fCs, fFlag, slice) \
       ((am)->funcs->map((am), (fsrName), (fStart), (fEnd), (fStrand), (fCs), (fFlag), (slice)))

#define ChainedAssemblyMapper_fastMap(am, fsrName, fStart, fEnd, fStrand, fCs, slice) \
       ((am)->funcs->fastMap((am), (fsrName), (fStart), (fEnd), (fStrand), (fCs), (slice)))

#define ChainedAssemblyMapper_listIds(am, fsrName, fStart, fEnd, fCs) \
       ((am)->funcs->listIds((am), (fsrName), (fStart), (fEnd), (fCs)))

#define ChainedAssemblyMapper_listSeqRegions(am, fsrName, fStart, fEnd, fCs) \
       ((am)->funcs->listSeqRegions((am), (fsrName), (fStart), (fEnd), (fCs)))

#define ChainedAssemblyMapper_mapCoordinatesToAssembly(am, contigName, fStart, fEnd, fStrand) \
       ((am)->funcs->mapCoordinatesToAssembly((am), (contigName), (fStart), (fEnd), (fStrand)))

#define ChainedAssemblyMapper_fastToAssembly(am, contigName, fStart, fEnd, fStrand) \
       ((am)->funcs->fastToAssembly((am), (contigName), (fStart), (fEnd), (fStrand)))

#define ChainedAssemblyMapper_mapCoordinatesToRawContig(am, chrName, fStart, fEnd, fStrand) \
       ((am)->funcs->mapCoordinatesToRawContig((am), (chrName), (fStart), (fEnd), (fStrand)))

#define ChainedAssemblyMapper_listContigIds(am, chrName, fStart, fEnd, fStrand) \
       ((am)->funcs->listContigIds((am), (chrName), (fStart), (fEnd), (fStrand)))

ChainedAssemblyMapper *ChainedAssemblyMapper_new(AssemblyMapperAdaptor *ama, Vector *coordSystems);

void ChainedAssemblyMapper_registerAllImpl(ChainedAssemblyMapper *cam);
void ChainedAssemblyMapper_flushImpl(ChainedAssemblyMapper *cam);
void ChainedAssemblyMapper_freeImpl(ChainedAssemblyMapper *cam);
int ChainedAssemblyMapper_getSize(ChainedAssemblyMapper *cam);

MapperRangeSet *ChainedAssemblyMapper_mapImpl(ChainedAssemblyMapper *cam, char *frmSeqRegionName, long frmStart, long frmEnd, int frmStrand,
                              CoordSystem *frmCs, int fastmap, Slice *toSlice);
MapperRangeSet *ChainedAssemblyMapper_fastMapImpl(ChainedAssemblyMapper *cam, char *frmSeqRegionName, long frmStart, long frmEnd, int frmStrand,
                              CoordSystem *frmCs);

MapperRangeSet *ChainedAssemblyMapper_mapCoordinatesToAssemblyImpl(ChainedAssemblyMapper *cam, char *contigName, long start, long end, int strand);
MapperRangeSet *ChainedAssemblyMapper_fastToAssemblyImpl(ChainedAssemblyMapper *cam, char *contigName, long start, long end, int strand);
MapperRangeSet *ChainedAssemblyMapper_mapCoordinatesToRawContigImpl(ChainedAssemblyMapper *cam, char *chrName, long start, long end, int strand);
Vector *ChainedAssemblyMapper_listContigIdsImpl(ChainedAssemblyMapper *cam, char *chrName, long start, long end, int strand);





Vector *ChainedAssemblyMapper_listIdsImpl(ChainedAssemblyMapper *cam, char *frmSeqRegionName, long frmStart, long frmEnd, CoordSystem *frmCs);
Vector *ChainedAssemblyMapper_listSeqRegionsImpl(ChainedAssemblyMapper *cam, char *frmSeqRegionName, long frmStart, long frmEnd, CoordSystem *frmCs);

#ifdef __CHAINEDASSEMBLYMAPPER_MAIN__
  ChainedAssemblyMapperFuncs
    chainedAssemblyMapperFuncs = {
                      ChainedAssemblyMapper_freeImpl, // free
                      ChainedAssemblyMapper_mapImpl,   // map
                      ChainedAssemblyMapper_flushImpl,   // flush
                      ChainedAssemblyMapper_fastMapImpl,  //fastMap
                      ChainedAssemblyMapper_listSeqRegionsImpl,  // listSeqRegions
                      ChainedAssemblyMapper_listIdsImpl,  // listIds
                      ChainedAssemblyMapper_registerAllImpl,  // registerAll
                      ChainedAssemblyMapper_mapCoordinatesToAssemblyImpl,  // mapCoordinatesToAssembly
                      ChainedAssemblyMapper_fastToAssemblyImpl,  // fastToAssembly
                      ChainedAssemblyMapper_mapCoordinatesToRawContigImpl,  // mapCoordinatesToRawContig
                      ChainedAssemblyMapper_listContigIdsImpl  // listContigIds
                   };
#else
  extern ChainedAssemblyMapperFuncs chainedAssemblyMapperFuncs;
#endif



#endif


