#ifndef __EXON_H__
#define __EXON_H__

#include <time.h>

#include "EnsC.h"
#include "DataModelTypes.h"
#include "AnnotatedSeqFeature.h"
#include "StableIdInfo.h"
#include "FeatureSet.h"
#include "Slice.h"

typedef struct ExonFuncsStruct {
  ANNOTATEDSEQFEATUREFUNCS_DATA
} ExonFuncs;

#define EXON_DATA \
  ANNOTATEDSEQFEATURE_DATA \
  FeatureSet components; \
  int stickyRank; \
  FeatureSet supportingFeatures;

#define FUNCSTRUCTTYPE ExonFuncs
struct ExonStruct {
  EXON_DATA
};
#undef FUNCSTRUCTTYPE

#define Exon_setStart(exon,start) AnnotatedSeqFeature_setStart((exon),start)
#define Exon_getStart(exon) AnnotatedSeqFeature_getStart((exon))

#define Exon_setEnd(exon,end) AnnotatedSeqFeature_setEnd((exon),end)
#define Exon_getEnd(exon) AnnotatedSeqFeature_getEnd((exon))

#define Exon_setScore(exon,score) AnnotatedSeqFeature_setScore((exon),score)
#define Exon_getScore(exon) AnnotatedSeqFeature_getScore((exon))

#define Exon_setpValue(exon,pValue) AnnotatedSeqFeature_setEValue((exon),pValue)
#define Exon_getpValue(exon) AnnotatedSeqFeature_getEValue((exon))

#define Exon_setPhase(exon,p) AnnotatedSeqFeature_setPhase((exon),(p))
#define Exon_getPhase(exon) AnnotatedSeqFeature_getPhase((exon))

#define Exon_setEndPhase(exon,ep) AnnotatedSeqFeature_setEndPhase((exon),(ep))
#define Exon_getEndPhase(exon) AnnotatedSeqFeature_getEndPhase((exon))

#define Exon_setStrand(exon,strand) AnnotatedSeqFeature_setStrand((exon),(strand))
#define Exon_getStrand(exon) AnnotatedSeqFeature_getStrand((exon))

#define Exon_setStableId(exon,stableId)  StableIdInfo_setStableId(&((exon)->si),(stableId))
char *Exon_getStableId(Exon *exon);

#define Exon_setVersion(exon,ver)  StableIdInfo_setVersion(&((exon)->si),(ver))
int Exon_getVersion(Exon *exon);

#define Exon_setCreated(exon,cd)  StableIdInfo_setCreated(&((exon)->si),(cd))
time_t Exon_getCreated(Exon *exon);

#define Exon_setModified(exon,mod)  StableIdInfo_setModified(&((exon)->si),(mod))
time_t Exon_getModified(Exon *exon);

#define Exon_setDbID(exon,dbID) AnnotatedSeqFeature_setDbID((exon),(dbID))
#define Exon_getDbID(exon) AnnotatedSeqFeature_getDbID((exon))

#define Exon_setAdaptor(exon,ad) AnnotatedSeqFeature_setAdaptor((exon),(ad))
#define Exon_getAdaptor(exon) AnnotatedSeqFeature_getAdaptor((exon))

#define Exon_setAnalysis(exon,ana) AnnotatedSeqFeature_setAnalysis((exon),(ana))
#define Exon_getAnalysis(exon) AnnotatedSeqFeature_getAnalysis((exon))

#define Exon_setStickyRank(exon,sr) (exon)->stickyRank = (sr)
#define Exon_getStickyRank(exon) (exon)->stickyRank

#define Exon_getLength(exon) AnnotatedSeqFeature_getLength((exon))

#define Exon_addComponentExon(exon, comp) FeatureSet_addFeature(&((exon)->components),(comp))

#define Exon_isSticky(exon) FeatureSet_getNumFeature(&((exon)->components))
#define Exon_getComponentExonAt(exon,ind) FeatureSet_getFeatureAt(&((exon)->components),(ind))
#define Exon_getComponents(exon) FeatureSet_getFeatures(&((exon)->components))
#define Exon_getNumComponentExon(exon) FeatureSet_getNumFeature(&((exon)->components))

#define Exon_addSupportingFeature(exon, sf) FeatureSet_addFeature(&((exon)->supportingFeatures),(sf))
#define Exon_getSupportingFeatureAt(exon,ind) FeatureSet_getFeatureAt(&((exon)->supportingFeatures),(ind))
#define Exon_getSupportingFeatures(exon) FeatureSet_getFeatures(&((exon)->supportingFeatures))
#define Exon_getNumSupportingFeature(exon) FeatureSet_getNumFeature(&((exon)->supportingFeatures))

void Exon_sortByStickyRank(Exon *exon);


#define Exon_setContig(exon,c) AnnotatedSeqFeature_setContig((exon),(c))
#define Exon_getContig(exon) AnnotatedSeqFeature_getContig((exon))

Exon *Exon_transformToSlice(Exon *exon, Slice *slice);


Exon *Exon_new();
Exon *Exon_copy(Exon *orig, CopyDepth depth);

#ifdef __EXON_MAIN__
  ExonFuncs exonFuncs = {NULL,NULL};
#else
  extern ExonFuncs exonFuncs;
#endif

#endif
