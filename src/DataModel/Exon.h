#ifndef __EXON_H__
#define __EXON_H__

#include <time.h>

#include "EnsC.h"
#include "DataModelTypes.h"
#include "SeqFeature.h"
#include "StableIdInfo.h"
#include "FeatureSet.h"
#include "Slice.h"

struct ExonStruct {
  SeqFeature sf;
  StableIdInfo si;
  FeatureSet components;
  int stickyRank;
};

#define Exon_setStart(exon,start) SeqFeature_setStart(&((exon)->sf),start)
#define Exon_getStart(exon) SeqFeature_getStart(&((exon)->sf))

#define Exon_setEnd(exon,end) SeqFeature_setEnd(&((exon)->sf),end)
#define Exon_getEnd(exon) SeqFeature_getEnd(&((exon)->sf))

#define Exon_setPhase(exon,p) SeqFeature_setPhase(&((exon)->sf),(p))
#define Exon_getPhase(exon) SeqFeature_getPhase(&((exon)->sf))

#define Exon_setEndPhase(exon,ep) SeqFeature_setEndPhase(&((exon)->sf),(ep))
#define Exon_getEndPhase(exon) SeqFeature_getEndPhase(&((exon)->sf))

#define Exon_setStrand(exon,strand) SeqFeature_setStrand(&((exon)->sf),(strand))
#define Exon_getStrand(exon) SeqFeature_getStrand(&((exon)->sf))

#define Exon_setStableId(exon,stableId)  StableIdInfo_setStableId(&((exon)->si),(stableId))
char *Exon_getStableId(Exon *exon);

#define Exon_setVersion(exon,ver)  StableIdInfo_setVersion(&((exon)->si),(ver))
int Exon_getVersion(Exon *exon);

#define Exon_setCreated(exon,cd)  StableIdInfo_setCreated(&((exon)->si),(cd))
time_t Exon_getCreated(Exon *exon);

#define Exon_setModified(exon,mod)  StableIdInfo_setModified(&((exon)->si),(mod))
time_t Exon_getModified(Exon *exon);

#define Exon_setDbID(e,dbID) SeqFeature_setDbID(&((e)->sf),(dbID))
#define Exon_getDbID(e) SeqFeature_getDbID(&((e)->sf))

#define Exon_setAdaptor(e,ad) SeqFeature_setAdaptor(&((e)->sf),(ad))
#define Exon_getAdaptor(e) SeqFeature_getAdaptor(&((e)->sf))

#define Exon_setAnalysis(e,ana) SeqFeature_setAnalysis(&((e)->sf),(ana))
#define Exon_getAnalysis(e) SeqFeature_getAnalysis(&((e)->sf))

#define Exon_setStickyRank(e,sr) (e)->stickyRank = (sr)
#define Exon_getStickyRank(e) (e)->stickyRank

#define Exon_getLength(e) SeqFeature_getLength(&((e)->sf))

#define Exon_addComponentExon(e, comp) FeatureSet_addFeature(&((e)->components),(comp))

#define Exon_isSticky(e) FeatureSet_getNumFeature(&((e)->components))
#define Exon_getComponentExonAt(e,ind) FeatureSet_getFeatureAt(&((e)->components),(ind))
#define Exon_getComponents(e) FeatureSet_getFeatures(&((e)->components))

#define Exon_getNumComponentExon(e) FeatureSet_getNumFeature(&((e)->components))

void Exon_sortByStickyRank(Exon *exon);


#define Exon_setContig(exon,c) SeqFeature_setContig(&((exon)->sf),(c))
#define Exon_getContig(exon) SeqFeature_getContig(&((exon)->sf))

Exon *Exon_transformToSlice(Exon *exon, Slice *slice);


Exon *Exon_new();
Exon *Exon_copy(Exon *orig, CopyDepth depth);


#endif
