#ifndef __STICKYEXON_H__
#define __STICKYEXON_H__

#include "Exon.h"
#include "DataModelTypes.h"
#include "Slice.h"
#include "EnsC.h"

EXONFUNC_TYPES(StickyExon)

typedef struct StickyExonFuncsStruct {
  EXONFUNCS_DATA(StickyExon)
} StickyExonFuncs;

#define STICKYEXON_DATA \
  EXON_DATA \
  FeatureSet components;

#define FUNCSTRUCTTYPE StickyExonFuncs
struct StickyExonStruct {
  STICKYEXON_DATA
};
#undef FUNCSTRUCTTYPE

#define StickyExon_setStart(exon,start) Exon_setStart((exon),(start))
#define StickyExon_getStart(exon) Exon_getStart((exon))

#define StickyExon_setSeqCacheString(exon,seq) Exon_setSeqCacheString((exon),(seq))
#define StickyExon_getSeqCacheString(exon) Exon_getSeqCacheString((exon))

#define StickyExon_setEnd(exon,end) Exon_setEnd((exon),(end))
#define StickyExon_getEnd(exon) Exon_getEnd((exon))

#define StickyExon_setScore(exon,score) Exon_setScore((exon),(score))
#define StickyExon_getScore(exon) Exon_getScore((exon))

#define StickyExon_setStableId(exon,sid) Exon_setStableId((exon),(sid))
#define StickyExon_getStableId(exon) Exon_getStableId((exon))

#define StickyExon_setVersion(exon,ver) Exon_setVersion((exon),(ver))
#define StickyExon_getVersion(exon) Exon_getVersion((exon))

#define StickyExon_setModified(exon,mod) Exon_setModified((exon),(mod))
#define StickyExon_getModified(exon) Exon_getModified((exon))

#define StickyExon_setCreated(exon,cd) Exon_setCreated((exon),(cd))
#define StickyExon_getCreated(exon) Exon_getCreated((exon))

#define StickyExon_setpValue(exon,pValue) Exon_setEValue((exon),(pValue))
#define StickyExon_getpValue(exon) Exon_getEValue((exon))

#define StickyExon_setPhase(exon,p) Exon_setPhase((exon),(p))
#define StickyExon_getPhase(exon) Exon_getPhase((exon))

#define StickyExon_setEndPhase(exon,ep) Exon_setEndPhase((exon),(ep))
#define StickyExon_getEndPhase(exon) Exon_getEndPhase((exon))

#define StickyExon_setStrand(exon,strand) Exon_setStrand((exon),(strand))
#define StickyExon_getStrand(exon) Exon_getStrand((exon))

#define Exon_setStableId(exon,stableId)  StableIdInfo_setStableId(&((exon)->si),(stableId))
char *Exon_getStableId(Exon *exon);

#define Exon_setVersion(exon,ver)  StableIdInfo_setVersion(&((exon)->si),(ver))
int Exon_getVersion(Exon *exon);

#define Exon_setCreated(exon,cd)  StableIdInfo_setCreated(&((exon)->si),(cd))
time_t Exon_getCreated(Exon *exon);

#define Exon_setModified(exon,mod)  StableIdInfo_setModified(&((exon)->si),(mod))
time_t Exon_getModified(Exon *exon);

#define StickyExon_setDbID(exon,dbID) Exon_setDbID((exon),(dbID))
#define StickyExon_getDbID(exon) Exon_getDbID((exon))

#define StickyExon_setAdaptor(exon,ad) Exon_setAdaptor((exon),(ad))
#define StickyExon_getAdaptor(exon) Exon_getAdaptor((exon))

#define StickyExon_setAnalysis(exon,ana) Exon_setAnalysis((exon),(ana))
#define StickyExon_getAnalysis(exon) Exon_getAnalysis((exon))

#define StickyExon_setContig(exon,c) Exon_setContig((exon),(c))
#define StickyExon_getContig(exon) Exon_getContig((exon))

// Can't have sticky exons in sticky exons
#define StickyExon_setStickyRank(exon,sr) ERROR;
#define StickyExon_getStickyRank(exon) ERROR;

int StickyExon_getLength(StickyExon *exon);

#define StickyExon_addComponentExon(exon, comp) FeatureSet_addFeature(&((exon)->components),(comp))

#define StickyExon_getComponentExonAt(exon,ind) FeatureSet_getFeatureAt(&((exon)->components),(ind))
#define StickyExon_getComponents(exon) FeatureSet_getFeatures(&((exon)->components))
#define StickyExon_getComponentExonCount(exon) FeatureSet_getNumFeature(&((exon)->components))

void StickyExon_addSupportingFeatures(StickyExon *exon, Vector *features);
Vector *StickyExon_getAllSupportingFeatures(StickyExon *stickyExon);

//#define StickyExon_getSupportingFeatureAt(exon,ind) Exon_getSupportingFeatureAt((exon),(ind))
//#define StickyExon_getSupportingFeatureCount(exon) Exon_getSupportingFeatureCount((exon))

void StickyExon_sortByStickyRank(StickyExon *exon);

int StickyExon_stickyRankCompFunc(const void *a, const void *b);

Exon *StickyExon_transformRawContigToSlice(StickyExon *exon, Slice *slice);
void StickyExon_loadGenomicMapper(StickyExon *stickyExon, Mapper *mapper, IDType id, int start);
Exon *StickyExon_adjustStartEnd(StickyExon *stickyExon, int startAdjust, int endAdjust);
char *StickyExon_getSeqString(StickyExon *stickyExon);
void StickyExon_addSupportingFeatures(StickyExon *stickyExon, Vector *features);

char *StickyExon_getPeptide(StickyExon *stickyExon, Transcript *trans);

StickyExon *StickyExon_new();

#ifdef __STICKYEXON_MAIN__
  StickyExonFuncs 
    stickyExonFuncs = {
                       NULL, // getStart
                       NULL, // setStart
                       NULL, // getEnd
                       NULL,  // setEnd
                       NULL, // getStrand
                       NULL, // setStrand
                       NULL, // getSeq
                       NULL, // setSeq
                       NULL, // getLength
                       NULL, // reverseComplement
                       (StickyExon_TransformToRawContigFunc)SeqFeature_transformToRawContigImpl,
                       (StickyExon_TransformToSliceFunc)SeqFeature_transformToSliceImpl,
                       StickyExon_transformRawContigToSlice,
                       NULL, // transformSliceToRawContig - SHOULD NOT HAVE ONE OF THESE FOR STICKIES BECAUSE THEY SHOULD NEVER BE ON SLICES
                       NULL, // transformSliceToSlice - OR ONE OF THESE
                       StickyExon_loadGenomicMapper,
                       StickyExon_adjustStartEnd,
                       StickyExon_getPeptide,
                       StickyExon_addSupportingFeatures,
                       StickyExon_getAllSupportingFeatures,
                       StickyExon_getSeqString
                      };
#else
  extern StickyExonFuncs stickyExonFuncs;
#endif

#endif
