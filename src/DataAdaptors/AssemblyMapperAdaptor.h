#ifndef __ASSEMBLYMAPPERADAPTOR_H__
#define __ASSEMBLYMAPPERADAPTOR_H__

#include "GenomicRange.h"
#include "AssemblyMapper.h"

#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "CoordSystem.h"

#include "StringHash.h"
#include "SeqRegionRange.h"


struct AssemblyMapperAdaptorStruct {
  BASEADAPTOR_DATA
  StringHash *asmMapperCache;
  StringHash *srNameCache;
  IDHash *srIdCache;
  IDHash *multSeqIdCache;
};

AssemblyMapperAdaptor *AssemblyMapperAdaptor_new(DBAdaptor *dba);

void AssemblyMapperAdaptor_cacheSeqIdsWithMultAssemblies(AssemblyMapperAdaptor *ama);
AssemblyMapper *AssemblyMapperAdaptor_fetchByCoordSystems(AssemblyMapperAdaptor *ama, CoordSystem *cs1, CoordSystem *cs2);
SeqRegionRange *AssemblyMapperAdaptor_addToRangeVector(Vector *ranges, IDType id, long start, long end, char *name);
void AssemblyMapperAdaptor_registerAssembled(AssemblyMapperAdaptor *ama, AssemblyMapper *assm, IDType asmSeqRegion, long asmStart, long asmEnd);
IDType AssemblyMapperAdaptor_seqRegionNameToId(AssemblyMapperAdaptor *ama, char *srName, IDType csId);
char *AssemblyMapperAdaptor_seqRegionIdToName(AssemblyMapperAdaptor *ama, IDType srId);
void AssemblyMapperAdaptor_addToSrCaches(AssemblyMapperAdaptor *ama, char *regionName, IDType regionId, IDType csId, long regionLength);
void AssemblyMapperAdaptor_registerComponent(AssemblyMapperAdaptor *ama, AssemblyMapper *asmMapper, IDType cmpSeqRegion);
void AssemblyMapperAdaptor_registerChained(AssemblyMapperAdaptor *ama, ChainedAssemblyMapper *casmMapper,
                                           char *from, IDType seqRegionId, Vector *ranges, Slice *toSlice );
void AssemblyMapperAdaptor_registerChainedSpecial(AssemblyMapperAdaptor *ama, ChainedAssemblyMapper *casmMapper,
                                                  char *from, IDType seqRegionId, Vector *ranges, Slice *toSlice);
void AssemblyMapperAdaptor_registerAll(AssemblyMapperAdaptor *ama, AssemblyMapper *asmMapper);
void AssemblyMapperAdaptor_registerAllChained(AssemblyMapperAdaptor *ama, ChainedAssemblyMapper *casmMapper);
void AssemblyMapperAdaptor_buildCombinedMapper(AssemblyMapperAdaptor *ama, Vector *ranges, Mapper *startMidMapper,
                                              Mapper *endMidMapper, Mapper *combinedMapper, char *startName);
Vector *AssemblyMapperAdaptor_seqRegionsToIds(AssemblyMapperAdaptor *ama, CoordSystem *coordSystem, Vector *seqRegions);
Vector *AssemblyMapperAdaptor_seqIdsToRegions(AssemblyMapperAdaptor *ama, Vector *seqRegionIds);
void AssemblyMapperAdaptor_deleteCache(AssemblyMapperAdaptor *ama);
void AssemblyMapperAdaptor_registerRegion(AssemblyMapperAdaptor *ama, AssemblyMapper *asmMapper, char *type,
                                          char *chrName, long start, long end);
void AssemblyMapperAdaptor_registerContig(AssemblyMapperAdaptor *ama, AssemblyMapper *asmMapper, char *type, IDType contigId);



//char *makeMappingPathKey(Vector *path);
char *makeMappingPathKey(Vector *path, char *key);

AssemblyMapper *AssemblyMapperAdaptor_fetchByType(AssemblyMapperAdaptor *ama, char *type);

#endif
