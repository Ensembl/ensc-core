#ifndef __ANNOTATEDSEQFEATURE_H__
#define __ANNOTATEDSEQFEATURE_H__

#include <time.h>

#include "EnsC.h"
#include "DataModelTypes.h"
#include "SeqFeature.h"
#include "StableIdInfo.h"

#define ANNOTATEDSEQFEATUREFUNCS_DATA
  SEQFEATUREFUNCS_DATA

typedef struct AnnotatedSeqFeatureFuncsStruct {
  ANNOTATEDSEQFEATUREFUNCS_DATA
} AnnotatedSeqFeatureFuncs;

#define ANNOTATEDSEQFEATURE_DATA \
  SEQFEATURE_DATA \
  StableIdInfo si;

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

#define AnnotatedSeqFeature_setEValue(asf,pValue) SeqFeature_setEValue((asf),pValue)
#define AnnotatedSeqFeature_getEValue(asf) SeqFeature_getEValue((asf))

#define AnnotatedSeqFeature_setPhase(asf,p) SeqFeature_setPhase((asf),(p))
#define AnnotatedSeqFeature_getPhase(asf) SeqFeature_getPhase((asf))

#define AnnotatedSeqFeature_setEndPhase(asf,ep) SeqFeature_setEndPhase((asf),(ep))
#define AnnotatedSeqFeature_getEndPhase(asf) SeqFeature_getEndPhase((asf))

#define AnnotatedSeqFeature_setStrand(asf,strand) SeqFeature_setStrand((asf),(strand))
#define AnnotatedSeqFeature_getStrand(asf) SeqFeature_getStrand((asf))

#define AnnotatedSeqFeature_setStableId(asf,stableId)  StableIdInfo_setStableId(&((asf)->si),(stableId))
char *AnnotatedSeqFeature_getStableId(AnnotatedSeqFeature *asf);

#define AnnotatedSeqFeature_setVersion(asf,ver)  StableIdInfo_setVersion(&((asf)->si),(ver))
int AnnotatedSeqFeature_getVersion(AnnotatedSeqFeature *asf);

#define AnnotatedSeqFeature_setCreated(asf,cd)  StableIdInfo_setCreated(&((asf)->si),(cd))
time_t AnnotatedSeqFeature_getCreated(AnnotatedSeqFeature *asf);

#define AnnotatedSeqFeature_setModified(asf,mod)  StableIdInfo_setModified(&((asf)->si),(mod))
time_t AnnotatedSeqFeature_getModified(AnnotatedSeqFeature *asf);

#define AnnotatedSeqFeature_setDbID(asf,dbID) SeqFeature_setDbID((asf),(dbID))
#define AnnotatedSeqFeature_getDbID(asf) SeqFeature_getDbID((asf))

#define AnnotatedSeqFeature_setAdaptor(asf,ad) SeqFeature_setAdaptor((asf),(ad))
#define AnnotatedSeqFeature_getAdaptor(asf) SeqFeature_getAdaptor((asf))

#define AnnotatedSeqFeature_setAnalysis(asf,ana) SeqFeature_setAnalysis((asf),(ana))
#define AnnotatedSeqFeature_getAnalysis(asf) SeqFeature_getAnalysis((asf))

#define AnnotatedSeqFeature_getLength(asf) SeqFeature_getLength((asf))

#define AnnotatedSeqFeature_setContig(asf,c) SeqFeature_setContig((asf),(c))
#define AnnotatedSeqFeature_getContig(asf) SeqFeature_getContig((asf))

#endif
