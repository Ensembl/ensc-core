#ifndef __REPEATFEATURE_H__
#define __REPEATFEATURE_H__

#include "DataModelTypes.h"
#include "SeqFeature.h"
#include "RepeatConsensus.h"

#define FUNCSTRUCTTYPE SeqFeatureFuncs
struct RepeatFeatureStruct {
  SEQFEATURE_DATA 
  int hitStart;
  int hitEnd;
  RepeatConsensus *repeatConsensus;
};
#undef FUNCSTRUCTTYPE

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

#define RepeatFeature_setContig(repeat,c) SeqFeature_setContig((repeat),(c))
#define RepeatFeature_getContig(repeat) SeqFeature_getContig((repeat))

#define RepeatFeature_setScore(repeat,score) SeqFeature_setScore((repeat),(score))
#define RepeatFeature_getScore(repeat) SeqFeature_getScore((repeat))

#define RepeatFeature_setHitStart(repeat,st) (repeat)->hitStart = (st)
#define RepeatFeature_getHitStart(repeat) (repeat)->hitStart

#define RepeatFeature_setHitEnd(repeat,end) (repeat)->hitEnd = (end)
#define RepeatFeature_getHitEnd(repeat) (repeat)->hitEnd

#define RepeatFeature_setConsensus(repeat,rc) (repeat)->repeatConsensus = (rc)
#define RepeatFeature_getConsensus(repeat) (repeat)->repeatConsensus

RepeatFeature *RepeatFeature_new();

#endif
