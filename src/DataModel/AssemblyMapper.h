#ifndef __ASSEMBLYMAPPER_H__
#define __ASSEMBLYMAPPER_H__

#include "DataModelTypes.h"
#include "EnsC.h"

#include "Vector.h"
#include "IDHash.h"
#include "Mapper.h"
#include "AdaptorTypes.h"
#include "MapperRangeSet.h"

#define MAXCHUNK 2000

struct AssemblyMapperStruct {
  AssemblyMapperAdaptor *adaptor;
  Mapper *mapper;
  IDHash *componentRegister;
  IDHash *assembledRegister;
  CoordSystem *assembledCoordSystem;
  CoordSystem *componentCoordSystem;
  int maxPairCount;
};

#define AssemblyMapper_setChrChunkHash(am,h) (am)->chrChunkHash = (h)
#define AssemblyMapper_getChrChunkHash(am) (am)->chrChunkHash

#define AssemblyMapper_setAdaptor(am, ad) (am)->adaptor = (ad)
#define AssemblyMapper_getAdaptor(am) (am)->adaptor

#define AssemblyMapper_setAssembledCoordSystem(am, cs) (am)->assembledCoordSystem = (cs)
#define AssemblyMapper_getAssembledCoordSystem(am) (am)->assembledCoordSystem

#define AssemblyMapper_setComponentCoordSystem(am, cs) (am)->componentCoordSystem = (cs)
#define AssemblyMapper_getComponentCoordSystem(am) (am)->componentCoordSystem

#define AssemblyMapper_setMaxPairCount(am, mp) (am)->maxPairCount = (mp)
#define AssemblyMapper_getMaxPairCount(am) (am)->maxPairCount

#define AssemblyMapper_setMapper(am, m) (am)->mapper = (m)
#define AssemblyMapper_getMapper(am) (am)->mapper

/*
char *AssemblyMapper_setType(AssemblyMapper *am, char *type);
#define AssemblyMapper_getType(am) (am)->type
*/

#define AssemblyMapper_setComponentRegister(am, cr) (am)->componentRegister = (cr)
#define AssemblyMapper_getComponentRegister(am) (am)->componentRegister

#define AssemblyMapper_setAssembledRegister(am, ar) (am)->assembledRegister = (ar)
#define AssemblyMapper_getAssembledRegister(am) (am)->assembledRegister

AssemblyMapper *AssemblyMapper_new(AssemblyMapperAdaptor *ama, Vector *coordSystems);

IDType AssemblyMapper_getSeqRegionId(AssemblyMapper *am, char *seqRegionName, CoordSystem *cs);
int AssemblyMapper_haveRegisteredComponent(AssemblyMapper *am, IDType cmpSeqRegionId);
int AssemblyMapper_haveRegisteredAssembled(AssemblyMapper *am, IDType asmSeqRegionId, int chunkId);
#endif
