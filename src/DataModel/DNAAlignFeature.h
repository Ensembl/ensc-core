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

#ifndef __DNAALIGNFEATURE_H__
#define __DNAALIGNFEATURE_H__

#include "DataModelTypes.h"
#include "BaseAlignFeature.h"

BASEALIGNFEATUREFUNC_TYPES(DNAAlignFeature)

typedef struct DNAAlignFeatureFuncsStruct {
  BASEALIGNFEATUREFUNCS_DATA(DNAAlignFeature)
} DNAAlignFeatureFuncs;
  
#define DNAALIGNFEATURE_DATA \
  BASEALIGNFEATURE_DATA \
  IDType pairDNAAlignFeatureId;

#ifdef COMMENTED_OUT
ECOSTRING extraData;
#endif

#define FUNCSTRUCTTYPE DNAAlignFeatureFuncs
struct DNAAlignFeatureStruct {
  DNAALIGNFEATURE_DATA
};
#undef FUNCSTRUCTTYPE

#define DNAAlignFeature_isStored(daf, db) Storable_isStored(&((daf)->st), (db))

DNAAlignFeature *DNAAlignFeature_new(void);

#define DNAAlignFeature_addFlag(fp, f) SeqFeature_addFlag((fp), (f))
#define DNAAlignFeature_getFlags(fp) SeqFeature_getFlags(fp)
#define DNAAlignFeature_removeFlag(fp, f) SeqFeature_removeFlag((fp), (f))

#define DNAAlignFeature_setCigarString(fp, ciggy) BaseAlignFeature_setCigarString((BaseAlignFeature *)(fp), (ciggy))
#define DNAAlignFeature_getCigarString(fp) BaseAlignFeature_getCigarString((fp))

#define DNAAlignFeature_setDbName(fp, dbName) BaseAlignFeature_setDbName((BaseAlignFeature *)(fp), (dbName))
#define DNAAlignFeature_getDbName(fp) BaseAlignFeature_getDbName((fp))

#define DNAAlignFeature_setDbDisplayName(fp, dbDisplayName) BaseAlignFeature_setDbDisplayName((BaseAlignFeature *)(fp), (dbDisplayName))
#define DNAAlignFeature_getDbDisplayName(fp) BaseAlignFeature_getDbDisplayName((fp))

#define DNAAlignFeature_sethCoverage(fp, hCoverage) BaseAlignFeature_sethCoverage((BaseAlignFeature *)(fp), (hCoverage))
#define DNAAlignFeature_gethCoverage(fp) BaseAlignFeature_gethCoverage((fp))

#define DNAAlignFeature_setExternalDbID(fp, externalDbId) BaseAlignFeature_setExternalDbID((BaseAlignFeature *)(fp), (externalDbId))
#define DNAAlignFeature_getExternalDbID(fp) BaseAlignFeature_getExternalDbID((fp))

#define DNAAlignFeature_setPairDNAAlignFeatureId(fp, dbId) (fp)->pairDNAAlignFeatureId = (dbId)
#define DNAAlignFeature_getPairDNAAlignFeatureId(fp) fp->pairDNAAlignFeatureId

ECOSTRING DNAAlignFeature_setExtraData(DNAAlignFeature *fp, char *extraData);
#define DNAAlignFeature_getExtraData(fp)  (fp)->extraData

#define DNAAlignFeature_getUngappedFeatures(fp) BaseAlignFeature_getUngappedFeatures((fp))

#define DNAAlignFeature_setHitSeqName(fp,stableId)  BaseAlignFeature_setHitSeqName((BaseAlignFeature *)(fp),(stableId))
#define DNAAlignFeature_getHitSeqName(fp)  BaseAlignFeature_getHitSeqName((fp))

#define DNAAlignFeature_setStart(fp,start) BaseAlignFeature_setStart((fp),(start))
#define DNAAlignFeature_getStart(fp) BaseAlignFeature_getStart((fp))

#define DNAAlignFeature_setHitStart(fp,start) BaseAlignFeature_setHitStart((fp),(start))
#define DNAAlignFeature_getHitStart(fp) BaseAlignFeature_getHitStart((fp))

#define DNAAlignFeature_setEnd(fp,end) BaseAlignFeature_setEnd((fp),(end))
#define DNAAlignFeature_getEnd(fp) BaseAlignFeature_getEnd((fp))

#define DNAAlignFeature_setHitEnd(fp,end) BaseAlignFeature_setHitEnd((fp),(end))
#define DNAAlignFeature_getHitEnd(fp) BaseAlignFeature_getHitEnd((fp))

#define DNAAlignFeature_setStrand(fp,strand) BaseAlignFeature_setStrand((fp),(strand))
#define DNAAlignFeature_getStrand(fp) BaseAlignFeature_getStrand((fp))

#define DNAAlignFeature_setHitStrand(fp,strand) BaseAlignFeature_setHitStrand((fp),(strand))
#define DNAAlignFeature_getHitStrand(fp) BaseAlignFeature_getHitStrand((fp))

#define DNAAlignFeature_setDbID(fp,dbID) BaseAlignFeature_setDbID((fp),(dbID))
#define DNAAlignFeature_getDbID(fp) BaseAlignFeature_getDbID((fp))

#define DNAAlignFeature_setAnalysis(fp,anal) BaseAlignFeature_setAnalysis((fp),(anal))
#define DNAAlignFeature_getAnalysis(fp) BaseAlignFeature_getAnalysis((fp))

#define DNAAlignFeature_setSlice(fp,contig) BaseAlignFeature_setSlice((fp),(contig))
#define DNAAlignFeature_getSlice(fp) BaseAlignFeature_getSlice((fp))

#define DNAAlignFeature_setScore(fp,score) BaseAlignFeature_setScore((fp),(score))
#define DNAAlignFeature_getScore(fp) BaseAlignFeature_getScore((fp))

#define DNAAlignFeature_setpValue(fp,ev) BaseAlignFeature_setpValue((fp),(ev))
#define DNAAlignFeature_getpValue(fp) BaseAlignFeature_getpValue((fp))

#define DNAAlignFeature_setPercId(fp,pid) BaseAlignFeature_setPercId((fp),(pid))
#define DNAAlignFeature_getPercId(fp) BaseAlignFeature_getPercId((fp))

#define DNAAlignFeature_setSeqName(fp,name) BaseAlignFeature_setSeqName((fp),(name))
#define DNAAlignFeature_getSeqName(fp) BaseAlignFeature_getSeqName((fp))

#define DNAAlignFeature_setSpecies(fp,sp) BaseAlignFeature_setSpecies((fp),(sp))
#define DNAAlignFeature_getSpecies(fp) BaseAlignFeature_getSpecies((fp))

#define DNAAlignFeature_setHitSpecies(fp,sp) BaseAlignFeature_setHitSpecies((fp),(sp))
#define DNAAlignFeature_getHitSpecies(fp) BaseAlignFeature_getHitSpecies((fp))

#define DNAAlignFeature_setAdaptor(fp,ad) BaseAlignFeature_setAdaptor((fp),(ad))
#define DNAAlignFeature_getAdaptor(fp) BaseAlignFeature_getAdaptor((fp))


int DNAAlignFeature_getHitUnit(void);
int DNAAlignFeature_getQueryUnit(void);

#define DNAAlignFeature_transformToRawContig(fp) BaseAlignFeature_transformToRawContig((fp))
#define DNAAlignFeature_transformToSlice(fp,slice) BaseAlignFeature_transformToSlice((fp),(slice))

#define DNAAlignFeature_getSeqRegionStart(daf) SeqFeature_getSeqRegionStart((SeqFeature *)(daf))
#define DNAAlignFeature_getSeqRegionEnd(daf) SeqFeature_getSeqRegionEnd((SeqFeature *)(daf))
#define DNAAlignFeature_getSeqRegionStrand(daf) SeqFeature_getSeqRegionStrand((SeqFeature *)(daf))

#define DNAAlignFeature_free(fp) BaseAlignFeature_free((fp))
// Why go down the inheritance tree here??
#define DNAAlignFeature_shallowCopy(fp) DNAAlignFeature_shallowCopyImpl((fp))

#define DNAAlignFeature_deepCopy(fp) DNAAlignFeature_deepCopyImpl((fp))

#define DNAAlignFeature_getLength(fp) BaseAlignFeature_getLength(fp)

void DNAAlignFeature_freeImpl(DNAAlignFeature *daf);
DNAAlignFeature *DNAAlignFeature_shallowCopyImpl(DNAAlignFeature *daf);
DNAAlignFeature *DNAAlignFeature_deepCopyImpl(DNAAlignFeature *daf);


#ifdef __DNAALIGNFEATURE_MAIN__
  DNAAlignFeatureFuncs 
    dnaAlignFeatureFuncs = {
                             DNAAlignFeature_freeImpl, // free
                             DNAAlignFeature_shallowCopyImpl, // shallowCopy
                             DNAAlignFeature_deepCopyImpl, // deepCopy
                             NULL, // getStart
                             NULL, // setStart
                             NULL, // getEnd
                             NULL, // setEnd
                             NULL, // getStrand
                             NULL, // setStrand
                             NULL, // getSeq
                             NULL, // setSeq
                             NULL, // getLength
                             (DNAAlignFeature_ReverseComplementFunc)BaseAlignFeature_reverseComplementImpl,
                             (DNAAlignFeature_TransformToRawContigFunc)SeqFeature_transformToRawContigImpl,
                             (DNAAlignFeature_TransformToSliceFunc)SeqFeature_transformToSliceImpl,
                             (DNAAlignFeature_TransformRawContigToSliceFunc)SeqFeature_transformRawContigToSliceImpl,
                             (DNAAlignFeature_TransformSliceToRawContigFunc)BaseAlignFeature_transformSliceToRawContigImpl,
                             NULL, // transformSliceToSlice
                             DNAAlignFeature_getHitUnit,
                             DNAAlignFeature_getQueryUnit
                            };
#else 
  extern DNAAlignFeatureFuncs dnaAlignFeatureFuncs;
#endif

#endif
