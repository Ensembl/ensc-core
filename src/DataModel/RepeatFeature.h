#ifndef __REPEATFEATURE_H__
#define __REPEATFEATURE_H__

#include "DataModelTypes.h"
#include "SeqFeature.h"
#include "RepeatConsensus.h"

struct RepeatFeatureStruct {
  SeqFeature sf;
  int hitStart;
  int hitEnd;
  RepeatConsensus *repeatConsensus;
};

#define RepeatFeature_setStart(repeat,start) SeqFeature_setStart(&((repeat)->sf),start)
#define RepeatFeature_getStart(repeat) SeqFeature_getStart(&((repeat)->sf))

#define RepeatFeature_setEnd(repeat,end) SeqFeature_setEnd(&((repeat)->sf),end)
#define RepeatFeature_getEnd(repeat) SeqFeature_getEnd(&((repeat)->sf))

#define RepeatFeature_setPhase(repeat,p) SeqFeature_setPhase(&((repeat)->sf),(p))
#define RepeatFeature_getPhase(repeat) SeqFeature_getPhase(&((repeat)->sf))

#define RepeatFeature_setEndPhase(repeat,ep) SeqFeature_setEndPhase(&((repeat)->sf),(ep))
#define RepeatFeature_getEndPhase(repeat) SeqFeature_getEndPhase(&((repeat)->sf))

#define RepeatFeature_setStrand(repeat,strand) SeqFeature_setStrand(&((repeat)->sf),(strand))
#define RepeatFeature_getStrand(repeat) SeqFeature_getStrand(&((repeat)->sf))

#define RepeatFeature_setDbID(repeat,dbID) SeqFeature_setDbID(&((repeat)->sf),(dbID))
#define RepeatFeature_getDbID(repeat) SeqFeature_getDbID(&((repeat)->sf))

#define RepeatFeature_setAdaptor(repeat,ad) SeqFeature_setAdaptor(&((repeat)->sf),(ad))
#define RepeatFeature_getAdaptor(repeat) SeqFeature_getAdaptor(&((repeat)->sf))

#define RepeatFeature_setAnalysis(repeat,ana) SeqFeature_setAnalysis(&((repeat)->sf),(ana))
#define RepeatFeature_getAnalysis(repeat) SeqFeature_getAnalysis(&((repeat)->sf))

#define RepeatFeature_setContig(repeat,c) SeqFeature_setContig(&((repeat)->sf),(c))
#define RepeatFeature_getContig(repeat) SeqFeature_getContig(&((repeat)->sf))

#define RepeatFeature_setScore(repeat,score) SeqFeature_setScore(&((repeat)->sf),(score))
#define RepeatFeature_getScore(repeat) SeqFeature_getScore(&((repeat)->sf))

#define RepeatFeature_setHitStart(repeat,st) (repeat)->hitStart = (st)
#define RepeatFeature_getHitStart(repeat) (repeat)->hitStart

#define RepeatFeature_setHitEnd(repeat,end) (repeat)->hitEnd = (end)
#define RepeatFeature_getHitEnd(repeat) (repeat)->hitEnd

#define RepeatFeature_setConsensus(repeat,rc) (repeat)->repeatConsensus = (rc)
#define RepeatFeature_getConsensus(repeat) (repeat)->repeatConsensus

RepeatFeature *RepeatFeature_new();

#endif
