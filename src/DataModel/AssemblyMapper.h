/*
 * Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

#ifndef __ASSEMBLYMAPPER_H__
#define __ASSEMBLYMAPPER_H__

#include "DataModelTypes.h"
#include "EnsC.h"
#include "BaseAssemblyMapper.h"

#include "Vector.h"
#include "IDHash.h"
#include "StringHash.h"
#include "Mapper.h"
#include "AdaptorTypes.h"
#include "MapperRangeSet.h"

// NIY: remove MAXCHUNK
#define MAXCHUNK 2000

BASEASSEMBLYMAPPERFUNC_TYPES(AssemblyMapper)

typedef struct AssemblyMapperFuncsStruct {
  BASEASSEMBLYMAPPERFUNCS_DATA(AssemblyMapper)
} AssemblyMapperFuncs;


#define FUNCSTRUCTTYPE AssemblyMapperFuncs
struct AssemblyMapperStruct {
  BASEASSEMBLYMAPPER_DATA
  Mapper *mapper;
  IDHash *componentRegister;
  IDHash *assembledRegister;
  CoordSystem *assembledCoordSystem;
  CoordSystem *componentCoordSystem;
};
#undef FUNCSTRUCTTYPE

#define AssemblyMapper_setAdaptor(am, ad) BaseAssemblyMapper_setAdaptor(am,ad)
#define AssemblyMapper_getAdaptor(am) BaseAssemblyMapper_getAdaptor(am)

#define AssemblyMapper_setAssembledCoordSystem(am, cs) (am)->assembledCoordSystem = (cs)
#define AssemblyMapper_getAssembledCoordSystem(am) (am)->assembledCoordSystem

#define AssemblyMapper_setComponentCoordSystem(am, cs) (am)->componentCoordSystem = (cs)
#define AssemblyMapper_getComponentCoordSystem(am) (am)->componentCoordSystem

#define AssemblyMapper_setMaxPairCount(am, mp) (am)->maxPairCount = (mp)
#define AssemblyMapper_getMaxPairCount(am) (am)->maxPairCount

#define AssemblyMapper_setMapper(am, m) (am)->mapper = (m)
#define AssemblyMapper_getMapper(am) (am)->mapper

#define AssemblyMapper_setComponentRegister(am, cr) (am)->componentRegister = (cr)
#define AssemblyMapper_getComponentRegister(am) (am)->componentRegister

#define AssemblyMapper_setAssembledRegister(am, ar) (am)->assembledRegister = (ar)
#define AssemblyMapper_getAssembledRegister(am) (am)->assembledRegister

AssemblyMapper *AssemblyMapper_new(AssemblyMapperAdaptor *ama, Vector *coordSystems);

//IDType AssemblyMapper_getSeqRegionId(AssemblyMapper *am, char *seqRegionName, CoordSystem *cs);
#define AssemblyMapper_getSeqRegionId(am, srName, cs) BaseAssemblyMapper_getSeqRegionId((BaseAssemblyMapper *)(am), (srName), (cs))

#define AssemblyMapper_registerAll(am) ((am)->funcs->registerAll((am)))

#define AssemblyMapper_flush(am) ((am)->funcs->flush((am)))

#define AssemblyMapper_free(am) ((am)->funcs->free((am)))

#define AssemblyMapper_map(am, fsrName, fStart, fEnd, fStrand, fCs, fFlag, slice) \
       ((am)->funcs->map((am), (fsrName), (fStart), (fEnd), (fStrand), (fCs), (fFlag), (slice)))

#define AssemblyMapper_fastMap(am, fsrName, fStart, fEnd, fStrand, fCs, slice) \
       ((am)->funcs->fastMap((am), (fsrName), (fStart), (fEnd), (fStrand), (fCs), (slice)))

#define AssemblyMapper_listIds(am, fsrName, fStart, fEnd, fCs) \
       ((am)->funcs->listIds((am), (fsrName), (fStart), (fEnd), (fCs)))

#define AssemblyMapper_listSeqRegions(am, fsrName, fStart, fEnd, fCs) \
       ((am)->funcs->listSeqRegions((am), (fsrName), (fStart), (fEnd), (fCs)))

#define AssemblyMapper_mapCoordinatesToAssembly(am, contigName, fStart, fEnd, fStrand) \
       ((am)->funcs->mapCoordinatesToAssembly((am), (contigName), (fStart), (fEnd), (fStrand)))

#define AssemblyMapper_fastToAssembly(am, contigName, fStart, fEnd, fStrand) \
       ((am)->funcs->fastToAssembly((am), (contigName), (fStart), (fEnd), (fStrand)))

#define AssemblyMapper_mapCoordinatesToRawContig(am, chrName, fStart, fEnd, fStrand) \
       ((am)->funcs->mapCoordinatesToRawContig((am), (chrName), (fStart), (fEnd), (fStrand)))

#define AssemblyMapper_listContigIds(am, chrName, fStart, fEnd, fStrand) \
       ((am)->funcs->listContigIds((am), (chrName), (fStart), (fEnd), (fStrand)))

int AssemblyMapper_haveRegisteredComponent(AssemblyMapper *am, IDType cmpSeqRegionId);
int AssemblyMapper_haveRegisteredAssembled(AssemblyMapper *am, IDType asmSeqRegionId, int chunkId);

void AssemblyMapper_freeImpl(AssemblyMapper *am);
void AssemblyMapper_flushImpl(AssemblyMapper *am);
MapperRangeSet *AssemblyMapper_mapImpl(AssemblyMapper *am, char *frmSeqRegionName, long frmStart, long frmEnd, int frmStrand, CoordSystem *frmCs, int fastMap, Slice *toSlice);
Vector *AssemblyMapper_listSeqRegionsImpl(AssemblyMapper *am, char *frmSeqRegionName, long frmStart, long frmEnd, CoordSystem *frmCs);
Vector *AssemblyMapper_listIdsImpl(AssemblyMapper *am, char *frmSeqRegionName, long frmStart, long frmEnd, CoordSystem *frmCs);
MapperRangeSet *AssemblyMapper_mapCoordinatesToAssemblyImpl(AssemblyMapper *am, char *contigName, long start, long end, int strand);
MapperRangeSet *AssemblyMapper_fastToAssemblyImpl(AssemblyMapper *am, char *contigName, long start, long end, int strand);
MapperRangeSet *AssemblyMapper_mapCoordinatesToRawContigImpl(AssemblyMapper *am, char *chrName, long start, long end, int strand);
Vector *AssemblyMapper_listContigIdsImpl(AssemblyMapper *am, char *chrName, long start, long end, int strand);
MapperRangeSet *AssemblyMapper_fastMapImpl(AssemblyMapper *am, char *frmSeqRegionName, long frmStart, long frmEnd, int frmStrand, CoordSystem *frmCs, Slice *toSlice);
void AssemblyMapper_registerAllImpl(AssemblyMapper *am);

int AssemblyMapper_getSize(AssemblyMapper *am);

void AssemblyMapper_registerComponent(AssemblyMapper *am, IDType cmpSeqRegionId);
void AssemblyMapper_registerAssembled(AssemblyMapper *am, IDType asmSeqRegionId, int chunkId);


#ifdef __ASSEMBLYMAPPER_MAIN__
  AssemblyMapperFuncs
    assemblyMapperFuncs = {
                      AssemblyMapper_freeImpl, // free
                      NULL, // shallowCopy
                      NULL, // deepCopy
                      AssemblyMapper_mapImpl,   // map
                      AssemblyMapper_flushImpl,   // flush
                      AssemblyMapper_fastMapImpl,  //fastMap
                      AssemblyMapper_listSeqRegionsImpl,  // listSeqRegions
                      AssemblyMapper_listIdsImpl,  // listIds
                      AssemblyMapper_registerAllImpl,  // registerAll
                      AssemblyMapper_mapCoordinatesToAssemblyImpl,  // mapCoordinatesToAssembly
                      AssemblyMapper_fastToAssemblyImpl,  // fastToAssembly
                      AssemblyMapper_mapCoordinatesToRawContigImpl,  // mapCoordinatesToRawContig
                      AssemblyMapper_listContigIdsImpl  // listContigIds
                   };
#else
  extern AssemblyMapperFuncs assemblyMapperFuncs;
#endif

#endif
