#ifndef __SIMPLEFEATURE_H__
#define __SIMPLEFEATURE_H__

#include "EnsC.h"
#include "DataModelTypes.h"
#include "SeqFeature.h"
#include "Slice.h"

struct SimpleFeatureStruct {
  SeqFeature sf;
  char *displayLabel;
};

#define SimpleFeature_setStart(simpleFeature,start) SeqFeature_setStart(&((simpleFeature)->sf),start)
#define SimpleFeature_getStart(simpleFeature) SeqFeature_getStart(&((simpleFeature)->sf))

#define SimpleFeature_setEnd(simpleFeature,end) SeqFeature_setEnd(&((simpleFeature)->sf),end)
#define SimpleFeature_getEnd(simpleFeature) SeqFeature_getEnd(&((simpleFeature)->sf))

#define SimpleFeature_setScore(simpleFeature,score) SeqFeature_setScore(&((simpleFeature)->sf),score)
#define SimpleFeature_getScore(simpleFeature) SeqFeature_getScore(&((simpleFeature)->sf))

#define SimpleFeature_setPhase(simpleFeature,p) SeqFeature_setPhase(&((simpleFeature)->sf),(p))
#define SimpleFeature_getPhase(simpleFeature) SeqFeature_getPhase(&((simpleFeature)->sf))

#define SimpleFeature_setEndPhase(simpleFeature,ep) SeqFeature_setEndPhase(&((simpleFeature)->sf),(ep))
#define SimpleFeature_getEndPhase(simpleFeature) SeqFeature_getEndPhase(&((simpleFeature)->sf))

#define SimpleFeature_setStrand(simpleFeature,strand) SeqFeature_setStrand(&((simpleFeature)->sf),(strand))
#define SimpleFeature_getStrand(simpleFeature) SeqFeature_getStrand(&((simpleFeature)->sf))

#define SimpleFeature_setDbID(simpleFeature,dbID) SeqFeature_setDbID(&((simpleFeature)->sf),(dbID))
#define SimpleFeature_getDbID(simpleFeature) SeqFeature_getDbID(&((simpleFeature)->sf))

#define SimpleFeature_setAdaptor(simpleFeature,ad) SeqFeature_setAdaptor(&((simpleFeature)->sf),(ad))
#define SimpleFeature_getAdaptor(simpleFeature) SeqFeature_getAdaptor(&((simpleFeature)->sf))

#define SimpleFeature_setAnalysis(simpleFeature,ana) SeqFeature_setAnalysis(&((simpleFeature)->sf),(ana))
#define SimpleFeature_getAnalysis(simpleFeature) SeqFeature_getAnalysis(&((simpleFeature)->sf))

#define SimpleFeature_getLength(simpleFeature) SeqFeature_getLength(&((simpleFeature)->sf))

#define SimpleFeature_setContig(simpleFeature,c) SeqFeature_setContig(&((simpleFeature)->sf),(c))
#define SimpleFeature_getContig(simpleFeature) SeqFeature_getContig(&((simpleFeature)->sf))

SimpleFeature *SimpleFeature_transformToSlice(SimpleFeature *simpleFeature, Slice *slice);

SimpleFeature *SimpleFeature_new();

char *SimpleFeature_setDisplayLabel(SimpleFeature *sf, char *label);
#define SimpleFeature_getDisplayLabel(simpleFeature) (simpleFeature)->displayLabel


#endif
