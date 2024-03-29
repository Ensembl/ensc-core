/*
 * Copyright [1999-2024] EMBL-European Bioinformatics Institute
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

#ifndef __ANNOTATEDSEQFEATURE_H__
#define __ANNOTATEDSEQFEATURE_H__

#include <time.h>

#include "EnsC.h"
#include "DataModelTypes.h"
#include "SeqFeature.h"
#include "StableIdInfo.h"
#include "DBEntry.h"

#define ANNOTATEDSEQFEATUREFUNC_TYPES(CLASSTYPE) \
  SEQFEATUREFUNC_TYPES(CLASSTYPE)

#define ANNOTATEDSEQFEATUREFUNCS_DATA(CLASSTYPE) \
  SEQFEATUREFUNCS_DATA(CLASSTYPE)

typedef struct AnnotatedSeqFeatureFuncsStruct {
  ANNOTATEDSEQFEATUREFUNCS_DATA(SeqFeature)
} AnnotatedSeqFeatureFuncs;

#define ANNOTATEDSEQFEATURE_DATA \
  SEQFEATURE_DATA \
  StableIdInfo si; \
  DBEntry *displayXref;

#define FUNCSTRUCTTYPE AnnotatedSeqFeatureFuncs
struct AnnotatedSeqFeatureStruct {
  ANNOTATEDSEQFEATURE_DATA
};
#undef FUNCSTRUCTTYPE

#define AnnotatedSeqFeature_setStart(asf,start) SeqFeature_setStart((asf),start)
#define AnnotatedSeqFeature_getStart(asf) SeqFeature_getStart((asf))

#define AnnotatedSeqFeature_setEnd(asf,end) SeqFeature_setEnd((asf),end)
#define AnnotatedSeqFeature_getEnd(asf) SeqFeature_getEnd((asf))

#define AnnotatedSeqFeature_setScore(asf,score) SeqFeature_setScore((asf),score)
#define AnnotatedSeqFeature_getScore(asf) SeqFeature_getScore((asf))

#define AnnotatedSeqFeature_setpValue(asf,pValue) SeqFeature_setpValue((asf),pValue)
#define AnnotatedSeqFeature_getpValue(asf) SeqFeature_getpValue((asf))

#define AnnotatedSeqFeature_setPhase(asf,p) SeqFeature_setPhase((asf),(p))
#define AnnotatedSeqFeature_getPhase(asf) SeqFeature_getPhase((asf))

#define AnnotatedSeqFeature_setEndPhase(asf,ep) SeqFeature_setEndPhase((asf),(ep))
#define AnnotatedSeqFeature_getEndPhase(asf) SeqFeature_getEndPhase((asf))

#define AnnotatedSeqFeature_setStrand(asf,strand) SeqFeature_setStrand((asf),(strand))
#define AnnotatedSeqFeature_getStrand(asf) SeqFeature_getStrand((asf))

#define AnnotatedSeqFeature_setStableId(asf,stableId)  StableIdInfo_setStableId(&((asf)->si),(stableId))
#define AnnotatedSeqFeature_getStableId(asf)  StableIdInfo_getStableId(&((asf)->si))

#define AnnotatedSeqFeature_setVersion(asf,ver)  StableIdInfo_setVersion(&((asf)->si),(ver))
int AnnotatedSeqFeature_getVersion(AnnotatedSeqFeature *asf);


#define AnnotatedSeqFeature_setCreated(asf,cd)  StableIdInfo_setCreated(&((asf)->si),(cd))
time_t AnnotatedSeqFeature_getCreated(AnnotatedSeqFeature *asf);

#define AnnotatedSeqFeature_setModified(asf,mod)  StableIdInfo_setModified(&((asf)->si),(mod))
time_t AnnotatedSeqFeature_getModified(AnnotatedSeqFeature *asf);

#define AnnotatedSeqFeature_setDisplayXref(asf,xref)  (asf)->displayXref = (xref)
#define AnnotatedSeqFeature_getDisplayXref(asf)  (asf)->displayXref
//DBEntry *AnnotatedSeqFeature_getDisplayXref(AnnotatedSeqFeature *asf);

#define AnnotatedSeqFeature_setDbID(asf,dbID) SeqFeature_setDbID((asf),(dbID))
#define AnnotatedSeqFeature_getDbID(asf) SeqFeature_getDbID((asf))

#define AnnotatedSeqFeature_setAdaptor(asf,ad) SeqFeature_setAdaptor((asf),(ad))
#define AnnotatedSeqFeature_getAdaptor(asf) SeqFeature_getAdaptor((asf))

#define AnnotatedSeqFeature_setAnalysis(asf,ana) SeqFeature_setAnalysis((asf),(ana))
#define AnnotatedSeqFeature_getAnalysis(asf) SeqFeature_getAnalysis((asf))

#define AnnotatedSeqFeature_getLength(asf) SeqFeature_getLength((asf))

#define AnnotatedSeqFeature_setContig(asf,c) SeqFeature_setContig((asf),(c))
#define AnnotatedSeqFeature_getContig(asf) SeqFeature_getContig((asf))

#define AnnotatedSeqFeature_setSlice(asf,c) SeqFeature_setSlice((asf),(c))
#define AnnotatedSeqFeature_getSlice(asf) SeqFeature_getSlice((asf))

#define AnnotatedSeqFeature_transformToSlice(asf, slice) SeqFeature_transformToSlice((asf),(slice))
#define AnnotatedSeqFeature_transformToRawContig(asf) SeqFeature_transformToRawContig((asf))

#define AnnotatedSeqFeature_free(asf) SeqFeature_free((asf))

#define AnnotatedSeqFeature_shallowCopy(asf) SeqFeature_shallowCopy((asf))

#ifdef __ANNOTATEDSEQFEATURE_MAIN__
  AnnotatedSeqFeatureFuncs 
       annotatedSeqFeatureFuncs = {
                                   NULL, // free
                                   NULL, // shallowCopy
                                   NULL, // deepCopy
                                   NULL, // getStart
                                   NULL, // setStart
                                   NULL, // getEnd
                                   NULL  // setEnd
                                  };
#else
  extern AnnotatedSeqFeatureFuncs annotatedSeqFeatureFuncs;
#endif

#endif
