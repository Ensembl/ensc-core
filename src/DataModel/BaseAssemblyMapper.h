/*
 * Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
 * Copyright [2016-2017] EMBL-European Bioinformatics Institute
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

#ifndef __BASEASSEMBLYMAPPER_H__
#define __BASEASSEMBLYMAPPER_H__


#include "DataModelTypes.h"
#include "EnsC.h"
#include "IDHash.h"
#include "StringHash.h"
#include "Storable.h"
#include "BaseContig.h"
#include "Vector.h"
#include "AdaptorTypes.h"
#include "MapperRangeSet.h"

#include "EnsRoot.h"

#define BASEASSEMBLYMAPPERFUNC_TYPES(CLASSTYPE) \
  OBJECTFUNC_TYPES(CLASSTYPE) \
  typedef MapperRangeSet * (*CLASSTYPE ## _MapFunc)(CLASSTYPE *am, char *frmSeqRegionName, long frmStart, long frmEnd, int frmStrand, CoordSystem *frmCs, int fakeFastMapFlag, Slice *toSlice); \
  typedef void (*CLASSTYPE ## _FlushFunc)(CLASSTYPE *am); \
  typedef MapperRangeSet *(*CLASSTYPE ## _FastMapFunc)(CLASSTYPE *am, char *frmSeqRegionName, long frmStart, long frmEnd, int frmStrand, CoordSystem *frmCs, Slice *toSlice); \
  typedef Vector * (*CLASSTYPE ## _ListSeqRegionsFunc)(CLASSTYPE *am, char *frmSeqRegionName, long frmStart, long frmEnd, CoordSystem *frmCs); \
  typedef Vector * (*CLASSTYPE ## _ListIdsFunc)(CLASSTYPE *am, char *frmSeqRegionName, long frmStart, long frmEnd, CoordSystem *frmCs); \
  typedef void (*CLASSTYPE ## _RegisterAllFunc)(CLASSTYPE *am); \
  typedef MapperRangeSet *(*CLASSTYPE ## _MapCoordinatesToAssemblyFunc)(CLASSTYPE *am, char *contigName, long start, long end, int strand); \
  typedef MapperRangeSet *(*CLASSTYPE ## _FastToAssemblyFunc)(CLASSTYPE *am, char *contigName, long start, long end, int strand); \
  typedef MapperRangeSet *(*CLASSTYPE ## _MapCoordinatesToRawContigFunc)(CLASSTYPE *am, char *chrName, long start, long end, int strand); \
  typedef Vector *(*CLASSTYPE ## _ListContigIdsFunc)(CLASSTYPE *am, char *chrName, long start, long end, int strand);

#define BASEASSEMBLYMAPPERFUNCS_DATA(CLASSTYPE) \
  OBJECTFUNCS_DATA(CLASSTYPE) \
  CLASSTYPE ## _MapFunc map; \
  CLASSTYPE ## _FlushFunc flush; \
  CLASSTYPE ## _FastMapFunc fastMap; \
  CLASSTYPE ## _ListSeqRegionsFunc listSeqRegions; \
  CLASSTYPE ## _ListIdsFunc listIds; \
  CLASSTYPE ## _RegisterAllFunc registerAll; \
  CLASSTYPE ## _MapCoordinatesToAssemblyFunc mapCoordinatesToAssembly; \
  CLASSTYPE ## _FastToAssemblyFunc fastToAssembly; \
  CLASSTYPE ## _MapCoordinatesToRawContigFunc mapCoordinatesToRawContig; \
  CLASSTYPE ## _ListContigIdsFunc listContigIds;

BASEASSEMBLYMAPPERFUNC_TYPES(BaseAssemblyMapper)

typedef struct BaseAssemblyMapperFuncsStruct {
  BASEASSEMBLYMAPPERFUNCS_DATA(BaseAssemblyMapper)
} BaseAssemblyMapperFuncs;

#define BASEASSEMBLYMAPPER_DATA \
  ENSROOT_DATA \
  AssemblyMapperAdaptor *adaptor; \
  int maxPairCount;

#define FUNCSTRUCTTYPE BaseAssemblyMapperFuncs
struct BaseAssemblyMapperStruct {
  BASEASSEMBLYMAPPER_DATA
};
#undef FUNCSTRUCTTYPE

IDType BaseAssemblyMapper_getSeqRegionId(BaseAssemblyMapper *am, char *seqRegionName, CoordSystem *cs);

#define BaseAssemblyMapper_setAdaptor(am, ad) (am)->adaptor = (ad)
#define BaseAssemblyMapper_getAdaptor(am) (am)->adaptor

#define BaseAssemblyMapper_registerAll(am) ((am)->funcs->registerAll((am)))

#define BaseAssemblyMapper_flush(am) ((am)->funcs->flush((am)))

#define BaseAssemblyMapper_free(am) ((am)->funcs->free((am)))

#define BaseAssemblyMapper_map(am, fsrName, fStart, fEnd, fStrand, fCs, fFlag, slice) \
       ((am)->funcs->map((am), (fsrName), (fStart), (fEnd), (fStrand), (fCs), (fFlag), (slice)))

#define BaseAssemblyMapper_fastMap(am, fsrName, fStart, fEnd, fStrand, fCs, slice) \
       ((am)->funcs->fastMap((am), (fsrName), (fStart), (fEnd), (fStrand), (fCs), (slice)))

#define BaseAssemblyMapper_listIds(am, fsrName, fStart, fEnd, fCs) \
       ((am)->funcs->listIds((am), (fsrName), (fStart), (fEnd), (fCs)))

#define BaseAssemblyMapper_listSeqRegions(am, fsrName, fStart, fEnd, fCs) \
       ((am)->funcs->listSeqRegions((am), (fsrName), (fStart), (fEnd), (fCs)))

#define BaseAssemblyMapper_mapCoordinatesToAssembly(am, contigName, fStart, fEnd, fStrand) \
       ((am)->funcs->mapCoordinatesToAssembly((am), (contigName), (fStart), (fEnd), (fStrand)))

#define BaseAssemblyMapper_fastToAssembly(am, contigName, fStart, fEnd, fStrand) \
       ((am)->funcs->fastToAssembly((am), (contigName), (fStart), (fEnd), (fStrand)))

#define BaseAssemblyMapper_mapCoordinatesToRawContig(am, chrName, fStart, fEnd, fStrand) \
       ((am)->funcs->mapCoordinatesToRawContig((am), (chrName), (fStart), (fEnd), (fStrand)))

#define BaseAssemblyMapper_listContigIds(am, chrName, fStart, fEnd, fStrand) \
       ((am)->funcs->listContigIds((am), (chrName), (fStart), (fEnd), (fStrand)))

void BaseAssemblyMapper_registerAllImpl(BaseAssemblyMapper *am);
void BaseAssemblyMapper_freeImpl(BaseAssemblyMapper *am);
void BaseAssemblyMapper_flushImpl(BaseAssemblyMapper *am);
MapperRangeSet *BaseAssemblyMapper_mapImpl(BaseAssemblyMapper *am, char *frmSeqRegionName, long frmStart, long frmEnd, int frmStrand,
                                           CoordSystem *frmCs, int fakeFastMapFlag, Slice *toSlice);
MapperRangeSet *BaseAssemblyMapper_fastMapImpl(BaseAssemblyMapper *am, char *frmSeqRegionName, long frmStart, long frmEnd, int frmStrand,
                                               CoordSystem *frmCs, Slice *toSlice);
Vector *BaseAssemblyMapper_listIdsImpl(BaseAssemblyMapper *am, char *frmSeqRegionName, long frmStart, long frmEnd, CoordSystem *frmCs);
Vector *BaseAssemblyMapper_listSeqRegionsImpl(BaseAssemblyMapper *am, char *frmSeqRegionName, long frmStart, long frmEnd, CoordSystem *frmCs);
MapperRangeSet *BaseAssemblyMapper_mapCoordinatesToAssemblyImpl(BaseAssemblyMapper *am, char *contigName, long start, long end, int strand);
MapperRangeSet *BaseAssemblyMapper_fastToAssemblyImpl(BaseAssemblyMapper *am, char *contigName, long start, long end, int strand);
MapperRangeSet *BaseAssemblyMapper_mapCoordinatesToRawContigImpl(BaseAssemblyMapper *am, char *chrName, long start, long end, int strand);
Vector *BaseAssemblyMapper_listContigIdsImpl(BaseAssemblyMapper *am, char *chrName, long start, long end, int strand);

#ifdef __BASEASSEMBLYMAPPER_MAIN__
 BaseAssemblyMapperFuncs 
   baseAssemblyMapperFuncs = {
                      BaseAssemblyMapper_freeImpl,   // free
                      NULL, // shallowCopy
                      NULL, // deepCopy
                      BaseAssemblyMapper_mapImpl,   // map
                      BaseAssemblyMapper_flushImpl,   // flush
                      BaseAssemblyMapper_fastMapImpl,  //fastMap
                      BaseAssemblyMapper_listSeqRegionsImpl,  // listSeqRegions
                      BaseAssemblyMapper_listIdsImpl,  // listIds
                      BaseAssemblyMapper_registerAllImpl,  // registerAll
                      BaseAssemblyMapper_mapCoordinatesToAssemblyImpl,  // mapCoordinatesToAssembly
                      BaseAssemblyMapper_fastToAssemblyImpl,  // fastToAssembly
                      BaseAssemblyMapper_mapCoordinatesToRawContigImpl,  // mapCoordinatesToRawContig
                      BaseAssemblyMapper_listContigIdsImpl  // listContigIds
                     };
#else
 extern BaseAssemblyMapperFuncs baseAssemblyMapperFuncs;
#endif

#endif
