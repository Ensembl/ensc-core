#ifndef __SIMPLEFEATURE_H__
#define __SIMPLEFEATURE_H__

#include "EnsC.h"
#include "DataModelTypes.h"
#include "SeqFeature.h"
#include "Slice.h"

#define SIMPLEFEATURE_DATA \
  SEQFEATURE_DATA \
  char *displayLabel;

#define FUNCSTRUCTTYPE SeqFeatureFuncs
struct SimpleFeatureStruct {
  SIMPLEFEATURE_DATA
};
#undef FUNCSTRUCTTYPE

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

#define SimpleFeature_setContig(simpleFeature,c) SeqFeature_setContig((simpleFeature),(c))
#define SimpleFeature_getContig(simpleFeature) SeqFeature_getContig((simpleFeature))

SimpleFeature *SimpleFeature_transformToSlice(SimpleFeature *simpleFeature, Slice *slice);

SimpleFeature *SimpleFeature_new();

char *SimpleFeature_setDisplayLabel(SimpleFeature *sf, char *label);
#define SimpleFeature_getDisplayLabel(simpleFeature) (simpleFeature)->displayLabel


#endif
