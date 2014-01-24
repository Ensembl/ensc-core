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

#ifndef __REPEATFEATURE_H__
#define __REPEATFEATURE_H__

#include "DataModelTypes.h"
#include "SeqFeature.h"
#include "RepeatConsensus.h"


SEQFEATUREFUNC_TYPES(RepeatFeature)

typedef struct RepeatFeatureFuncsStruct {
  SEQFEATUREFUNCS_DATA(RepeatFeature)
} RepeatFeatureFuncs;



#define FUNCSTRUCTTYPE RepeatFeatureFuncs
struct RepeatFeatureStruct {
  SEQFEATURE_DATA 
  int hitStart;
  int hitEnd;
  RepeatConsensus *repeatConsensus;
};
#undef FUNCSTRUCTTYPE

#define RepeatFeature_isStored(repeat, db) Storable_isStored(&((repeat)->st), (db))

#define RepeatFeature_getSeqRegionStart(repeat) SeqFeature_getSeqRegionStart((repeat))
#define RepeatFeature_getSeqRegionEnd(repeat) SeqFeature_getSeqRegionEnd((repeat))
#define RepeatFeature_getSeqRegionStrand(repeat) SeqFeature_getSeqRegionStrand((repeat))

#define RepeatFeature_setStart(repeat,start) SeqFeature_setStart((repeat),start)
#define RepeatFeature_getStart(repeat) SeqFeature_getStart((repeat))

#define RepeatFeature_setEnd(repeat,end) SeqFeature_setEnd((repeat),end)
#define RepeatFeature_getEnd(repeat) SeqFeature_getEnd((repeat))

#define RepeatFeature_setPhase(repeat,p) SeqFeature_setPhase((repeat),(p))
#define RepeatFeature_getPhase(repeat) SeqFeature_getPhase((repeat))

#define RepeatFeature_setEndPhase(repeat,ep) SeqFeature_setEndPhase((repeat),(ep))
#define RepeatFeature_getEndPhase(repeat) SeqFeature_getEndPhase((repeat))

#define RepeatFeature_setStrand(repeat,strand) SeqFeature_setStrand((repeat),(strand))
#define RepeatFeature_getStrand(repeat) SeqFeature_getStrand((repeat))

#define RepeatFeature_setDbID(repeat,dbID) SeqFeature_setDbID((repeat),(dbID))
#define RepeatFeature_getDbID(repeat) SeqFeature_getDbID((repeat))

#define RepeatFeature_setAdaptor(repeat,ad) SeqFeature_setAdaptor((repeat),(ad))
#define RepeatFeature_getAdaptor(repeat) SeqFeature_getAdaptor((repeat))

#define RepeatFeature_setAnalysis(repeat,ana) SeqFeature_setAnalysis((repeat),(ana))
#define RepeatFeature_getAnalysis(repeat) SeqFeature_getAnalysis((repeat))

#define RepeatFeature_setSlice(repeat,s) SeqFeature_setSlice((repeat),(s))
#define RepeatFeature_getSlice(repeat) SeqFeature_getSlice((repeat))

#define RepeatFeature_setScore(repeat,score) SeqFeature_setScore((repeat),(score))
#define RepeatFeature_getScore(repeat) SeqFeature_getScore((repeat))

#define RepeatFeature_setHitStart(repeat,st) (repeat)->hitStart = (st)
#define RepeatFeature_getHitStart(repeat) (repeat)->hitStart

#define RepeatFeature_setHitEnd(repeat,end) (repeat)->hitEnd = (end)
#define RepeatFeature_getHitEnd(repeat) (repeat)->hitEnd

#define RepeatFeature_setConsensus(repeat,rc) (repeat)->repeatConsensus = (rc)
#define RepeatFeature_getConsensus(repeat) (repeat)->repeatConsensus

#define RepeatFeature_transformToSlice(repeat,slice) SeqFeature_transformToSlice((repeat),(slice))
#define RepeatFeature_transformToRawContig(repeat) SeqFeature_transformToRawContig((repeat))

#define RepeatFeature_free(repeat) SeqFeature_free((repeat))

RepeatFeature *RepeatFeature_new();

void RepeatFeature_freeImpl(RepeatFeature *repeat);

#ifdef __REPEATFEATURE_MAIN__
 RepeatFeatureFuncs 
   repeatFeatureFuncs = {
                      RepeatFeature_freeImpl, // free
                      NULL, // shallowCopy
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
                      (RepeatFeature_TransformToRawContigFunc)SeqFeature_transformToRawContigImpl,
                      (RepeatFeature_TransformToSliceFunc)SeqFeature_transformToSliceImpl,
                      (RepeatFeature_TransformRawContigToSliceFunc)SeqFeature_transformRawContigToSliceImpl,
                      (RepeatFeature_TransformSliceToRawContigFunc)SeqFeature_transformSliceToRawContigImpl,
                      NULL // transformSliceToSlice
                     };
#else
 extern RepeatFeatureFuncs repeatFeatureFuncs;
#endif

#endif
