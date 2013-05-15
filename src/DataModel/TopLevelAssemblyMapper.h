#ifndef __TOPLEVELASSEMBLYMAPPER_H__
#define __TOPLEVELASSEMBLYMAPPER_H__

#include "BaseAssemblyMapper.h"
#include "DataModelTypes.h"
#include "EnsC.h"

#include "Vector.h"
#include "IDHash.h"
#include "Mapper.h"
#include "AdaptorTypes.h"
#include "MapperRangeSet.h"

#define MAXCHUNK 2000

BASEASSEMBLYMAPPERFUNC_TYPES(TopLevelAssemblyMapper)

typedef struct TopLevelAssemblyMapperFuncsStruct {
  BASEASSEMBLYMAPPERFUNCS_DATA(TopLevelAssemblyMapper)
} TopLevelAssemblyMapperFuncs;


#define FUNCSTRUCTTYPE TopLevelAssemblyMapperFuncs
struct TopLevelAssemblyMapperStruct {
  BASEASSEMBLYMAPPER_DATA
  Vector *coordSystems;
  CoordSystem *topLevelCs;
  CoordSystem *otherCs;
};
#undef FUNCSTRUCTTYPE

#define TopLevelAssemblyMapper_setAdaptor(tlam, ad) BaseAssemblyMapper_setAdaptor(tlam,ad)
#define TopLevelAssemblyMapper_getAdaptor(tlam) BaseAssemblyMapper_getAdaptor(tlam)

#define TopLevelAssemblyMapper_setTopLevelCoordSystem(tlam, tlcs) (tlam)->topLevelCs = tlcs
#define TopLevelAssemblyMapper_getTopLevelCoordSystem(tlam) (tlam)->topLevelCs

#define TopLevelAssemblyMapper_setOtherCoordSystem(tlam, ocs) (tlam)->otherCs = ocs
#define TopLevelAssemblyMapper_getOtherCoordSystem(tlam) (tlam)->otherCs

#define TopLevelAssemblyMapper_getCoordSystems(tlam) (tlam)->coordSystems

#define TopLevelAssemblyMapper_getSeqRegionId(am, srName, cs) BaseAssemblyMapper_getSeqRegionId((BaseAssemblyMapper *)(am), (srName), (cs))
//IDType TopLevelAssemblyMapper_getSeqRegionId(TopLevelAssemblyMapper *tlam, char *seqRegionName, CoordSystem *cs);

#define TopLevelAssemblyMapper_registerAll(am) BaseAssemblyMapper_registerAll((am))

#define TopLevelAssemblyMapper_flush(am) ((am)->funcs->flush((am)))

#define TopLevelAssemblyMapper_free(am) ((am)->funcs->free((am)))

#define TopLevelAssemblyMapper_map(am, fsrName, fStart, fEnd, fStrand, fCs, fFlag, slice) \
       ((am)->funcs->map((am), (fsrName), (fStart), (fEnd), (fStrand), (fCs), (fFlag), (slice)))

#define TopLevelAssemblyMapper_fastMap(am, fsrName, fStart, fEnd, fStrand, fCs, slice) \
       ((am)->funcs->fastMap((am), (fsrName), (fStart), (fEnd), (fStrand), (fCs), (slice)))

#define TopLevelAssemblyMapper_listIds(am, fsrName, fStart, fEnd, fCs) \
       ((am)->funcs->listIds((am), (fsrName), (fStart), (fEnd), (fCs)))

#define TopLevelAssemblyMapper_listSeqRegions(am, fsrName, fStart, fEnd, fCs) \
       ((am)->funcs->listSeqRegions((am), (fsrName), (fStart), (fEnd), (fCs)))

#define TopLevelAssemblyMapper_mapCoordinatesToAssembly(am, contigName, fStart, fEnd, fStrand) \
          BaseAssemblyMapper_mapCoordinatesToAssembly((am), (contigName), (fStart), (fEnd), (fStrand))

#define TopLevelAssemblyMapper_fastToAssembly(am, contigName, fStart, fEnd, fStrand) \
          BaseAssemblyMapper_fastToAssembly((am), (contigName), (fStart), (fEnd), (fStrand))

#define TopLevelAssemblyMapper_mapCoordinatesToRawContig(am, chrName, fStart, fEnd, fStrand) \
          BaseAssemblyMapper_mapCoordinatesToRawContig((am), (chrName), (fStart), (fEnd), (fStrand))

#define TopLevelAssemblyMapper_listContigIds(am, chrName, fStart, fEnd, fStrand) \
          BaseAssemblyMapper_listContigIds((am), (chrName), (fStart), (fEnd), (fStrand))

TopLevelAssemblyMapper *TopLevelAssemblyMapper_new(AssemblyMapperAdaptor *ama, CoordSystem *topLevelCs, CoordSystem *otherCs);

MapperRangeSet *TopLevelAssemblyMapper_mapImpl(TopLevelAssemblyMapper *tlam, char *frmSeqRegionName, long frmStart, long frmEnd, int frmStrand,
                                           CoordSystem *frmCs, int fastMap, Slice *fakeSliceArg);

void TopLevelAssemblyMapper_flushImpl(TopLevelAssemblyMapper *tlam);
MapperRangeSet *TopLevelAssemblyMapper_fastMapImpl(TopLevelAssemblyMapper *tlam, char *frmSeqRegionName, long frmStart, long frmEnd, int frmStrand,
                                               CoordSystem *frmCs, Slice *fakeSliceArg);
CoordSystem *TopLevelAssemblyMapper_assembledCoordSystem(TopLevelAssemblyMapper *tlam);
CoordSystem *TopLevelAssemblyMapper_componentCoordSystem(TopLevelAssemblyMapper *tlam);
Vector * TopLevelAssemblyMapper_list(TopLevelAssemblyMapper *tlam, char *frmSeqRegionName, long frmStart, long frmEnd, CoordSystem *frmCs, int SeqRegionsFlag);
Vector *TopLevelAssemblyMapper_listSeqRegionsImpl(TopLevelAssemblyMapper *tlam, char *frmSeqRegionName, long frmStart, long frmEnd, CoordSystem *frmCs);
Vector *TopLevelAssemblyMapper_listIdsImpl(TopLevelAssemblyMapper *tlam, char *frmSeqRegionName, long frmStart, long frmEnd, CoordSystem *frmCs);
void TopLevelAssemblyMapper_freeImpl(TopLevelAssemblyMapper *tlam);

void TopLevelAssemblyMapper_registerAllImpl(TopLevelAssemblyMapper *am);
MapperRangeSet *TopLevelAssemblyMapper_mapCoordinatesToAssemblyImpl(TopLevelAssemblyMapper *am, char *contigName, long start, long end, int strand);
MapperRangeSet *TopLevelAssemblyMapper_fastToAssemblyImpl(TopLevelAssemblyMapper *am, char *contigName, long start, long end, int strand);
MapperRangeSet *TopLevelAssemblyMapper_mapCoordinatesToRawContigImpl(TopLevelAssemblyMapper *am, char *chrName, long start, long end, int strand);
Vector *TopLevelAssemblyMapper_listContigIdsImpl(TopLevelAssemblyMapper *am, char *chrName, long start, long end, int strand);


#ifdef __TOPLEVELASSEMBLYMAPPER_MAIN__
  TopLevelAssemblyMapperFuncs
    topLevelAssemblyMapperFuncs = {
                      TopLevelAssemblyMapper_freeImpl, // free
                      NULL, // shallowCopy
                      NULL, // deepCopy
                      TopLevelAssemblyMapper_mapImpl,   // map
                      TopLevelAssemblyMapper_flushImpl,   // flush
                      TopLevelAssemblyMapper_fastMapImpl,  //fastMap
                      TopLevelAssemblyMapper_listSeqRegionsImpl,  // listSeqRegions
                      TopLevelAssemblyMapper_listIdsImpl,  // listIds
                      TopLevelAssemblyMapper_registerAllImpl,  // registerAll
                      TopLevelAssemblyMapper_mapCoordinatesToAssemblyImpl,  // mapCoordinatesToAssembly
                      TopLevelAssemblyMapper_fastToAssemblyImpl,  // fastToAssembly
                      TopLevelAssemblyMapper_mapCoordinatesToRawContigImpl,  // mapCoordinatesToRawContig
                      TopLevelAssemblyMapper_listContigIdsImpl  // listContigIds
                   };
#else
  extern TopLevelAssemblyMapperFuncs topLevelAssemblyMapperFuncs;
#endif


#endif
