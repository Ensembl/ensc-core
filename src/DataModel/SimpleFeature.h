/*
 * Copyright [1999-2017] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

#ifndef __SIMPLEFEATURE_H__
#define __SIMPLEFEATURE_H__

#include "EnsC.h"
#include "DataModelTypes.h"
#include "SeqFeature.h"
#include "Slice.h"

SEQFEATUREFUNC_TYPES(SimpleFeature)

typedef struct SimpleFeatureFuncsStruct {
  SEQFEATUREFUNCS_DATA(SimpleFeature)
} SimpleFeatureFuncs;



#define SIMPLEFEATURE_DATA \
  SEQFEATURE_DATA \
  char *displayLabel;

#define FUNCSTRUCTTYPE SimpleFeatureFuncs
struct SimpleFeatureStruct {
  SIMPLEFEATURE_DATA
};
#undef FUNCSTRUCTTYPE

#define SimpleFeature_isStored(simpleFeature, db) Storable_isStored(&((simpleFeature)->st), (db))

#define SimpleFeature_getSeqRegionStart(simpleFeature) SeqFeature_getSeqRegionStart((simpleFeature))
#define SimpleFeature_getSeqRegionEnd(simpleFeature) SeqFeature_getSeqRegionEnd((simpleFeature))
#define SimpleFeature_getSeqRegionStrand(simpleFeature) SeqFeature_getSeqRegionStrand((simpleFeature))

#define SimpleFeature_setStart(simpleFeature,start) SeqFeature_setStart((simpleFeature),start)
#define SimpleFeature_getStart(simpleFeature) SeqFeature_getStart((simpleFeature))

#define SimpleFeature_setEnd(simpleFeature,end) SeqFeature_setEnd((simpleFeature),end)
#define SimpleFeature_getEnd(simpleFeature) SeqFeature_getEnd((simpleFeature))

#define SimpleFeature_setScore(simpleFeature,score) SeqFeature_setScore((simpleFeature),score)
#define SimpleFeature_getScore(simpleFeature) SeqFeature_getScore((simpleFeature))

#define SimpleFeature_setPhase(simpleFeature,p) SeqFeature_setPhase((simpleFeature),(p))
#define SimpleFeature_getPhase(simpleFeature) SeqFeature_getPhase((simpleFeature))

#define SimpleFeature_setEndPhase(simpleFeature,ep) SeqFeature_setEndPhase((simpleFeature),(ep))
#define SimpleFeature_getEndPhase(simpleFeature) SeqFeature_getEndPhase((simpleFeature))

#define SimpleFeature_setStrand(simpleFeature,strand) SeqFeature_setStrand((simpleFeature),(strand))
#define SimpleFeature_getStrand(simpleFeature) SeqFeature_getStrand((simpleFeature))

#define SimpleFeature_setDbID(simpleFeature,dbID) SeqFeature_setDbID((simpleFeature),(dbID))
#define SimpleFeature_getDbID(simpleFeature) SeqFeature_getDbID((simpleFeature))

#define SimpleFeature_setAdaptor(simpleFeature,ad) SeqFeature_setAdaptor((simpleFeature),(ad))
#define SimpleFeature_getAdaptor(simpleFeature) SeqFeature_getAdaptor((simpleFeature))

#define SimpleFeature_setAnalysis(simpleFeature,ana) SeqFeature_setAnalysis((simpleFeature),(ana))
#define SimpleFeature_getAnalysis(simpleFeature) SeqFeature_getAnalysis((simpleFeature))

#define SimpleFeature_getLength(simpleFeature) SeqFeature_getLength((simpleFeature))

#define SimpleFeature_setSlice(simpleFeature,s) SeqFeature_setSlice((simpleFeature),(s))
#define SimpleFeature_getSlice(simpleFeature) SeqFeature_getSlice((simpleFeature))

#define SimpleFeature_transformToSlice(simpleFeature,slice) SeqFeature_transformToSlice((simpleFeature),(slice))
#define SimpleFeature_transformToRawContig(simpleFeature) SeqFeature_transformToRawContig((simpleFeature))

#define SimpleFeature_free(simpleFeature) SeqFeature_free((simpleFeature))
#define SimpleFeature_shallowCopy(simpleFeature) SeqFeature_shallowCopy((simpleFeature))

SimpleFeature *SimpleFeature_new();

ECOSTRING SimpleFeature_setDisplayLabel(SimpleFeature *sf, char *label);
#define SimpleFeature_getDisplayLabel(simpleFeature) (simpleFeature)->displayLabel

SimpleFeature *SimpleFeature_shallowCopyImpl(SimpleFeature *sf);

void SimpleFeature_freeImpl(SimpleFeature *sf);

#ifdef __SIMPLEFEATURE_MAIN__
 SimpleFeatureFuncs
   simpleFeatureFuncs = {
                      SimpleFeature_freeImpl,
                      SimpleFeature_shallowCopyImpl, // shallowCopy
                      NULL, // deepCopy
                      NULL, // getStart
                      NULL, // setStart
                      NULL, // getEnd
                      NULL, // setEnd
                      NULL, // getStrand
                      NULL, // setStrand
                      NULL, // getSeq
                      NULL, // setSeq
                      NULL, // getLength
                      NULL, // reverseComplement
                      (SimpleFeature_TransformToRawContigFunc)SeqFeature_transformToRawContigImpl,
                      (SimpleFeature_TransformToSliceFunc)SeqFeature_transformToSliceImpl,
                      (SimpleFeature_TransformRawContigToSliceFunc)SeqFeature_transformRawContigToSliceImpl,
                      (SimpleFeature_TransformSliceToRawContigFunc)SeqFeature_transformSliceToRawContigImpl,
                      NULL // transformSliceToSlice
                     };
#else
 extern SimpleFeatureFuncs simpleFeatureFuncs;
#endif


#endif
