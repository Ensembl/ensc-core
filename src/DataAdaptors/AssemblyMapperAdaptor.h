/*
 * Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
IDType AssemblyMapperAdaptor_seqRegionToId(AssemblyMapperAdaptor *ama, CoordSystem *coordSystem, char *seqRegionName);
Vector *AssemblyMapperAdaptor_seqIdsToRegions(AssemblyMapperAdaptor *ama, Vector *seqRegionIds);
void AssemblyMapperAdaptor_deleteCache(AssemblyMapperAdaptor *ama);
void AssemblyMapperAdaptor_registerRegion(AssemblyMapperAdaptor *ama, AssemblyMapper *asmMapper, char *type,
                                          char *chrName, long start, long end);
void AssemblyMapperAdaptor_registerContig(AssemblyMapperAdaptor *ama, AssemblyMapper *asmMapper, char *type, IDType contigId);



//char *makeMappingPathKey(Vector *path);
char *makeMappingPathKey(Vector *path, char *key);

AssemblyMapper *AssemblyMapperAdaptor_fetchByType(AssemblyMapperAdaptor *ama, char *type);

#endif
