#ifndef __EXON_H__
#define __EXON_H__

#include <time.h>

#include "EnsC.h"
#include "DataModelTypes.h"
#include "AnnotatedSeqFeature.h"
#include "StableIdInfo.h"
#include "FeatureSet.h"
#include "Slice.h"

#define EXONFUNC_TYPES(CLASSTYPE) \
  ANNOTATEDSEQFEATUREFUNC_TYPES(CLASSTYPE) \
  typedef void     (* CLASSTYPE ## _LoadGenomicMapperFunc)(CLASSTYPE *exon, Mapper *mapper, IDType id, int start); \
  typedef CLASSTYPE *   (* CLASSTYPE ## _AdjustStartEndFunc)(CLASSTYPE *exon, int startAdjust, int endAdjust); \
  typedef char *   (* CLASSTYPE ## _GetPeptideFunc)(CLASSTYPE *exon); \
  typedef void     (* CLASSTYPE ## _AddSupportingFeatureFunc)(CLASSTYPE *exon, SeqFeature *sf); \
  typedef Vector * (* CLASSTYPE ## _GetAllSupportingFeaturesFunc)(CLASSTYPE *exon);

EXONFUNC_TYPES(Exon)

#define EXONFUNCS_DATA(CLASSTYPE) \
  ANNOTATEDSEQFEATUREFUNCS_DATA(CLASSTYPE) \
  CLASSTYPE ## _LoadGenomicMapperFunc loadGenomicMapper; \
  CLASSTYPE ## _AdjustStartEndFunc adjustStartEnd; \
  CLASSTYPE ## _GetPeptideFunc getPeptide; \
  CLASSTYPE ## _AddSupportingFeatureFunc addSupportingFeature; \
  CLASSTYPE ## _GetAllSupportingFeaturesFunc getAllSupportingFeatures;

typedef struct ExonFuncsStruct {
  EXONFUNCS_DATA(Exon)
} ExonFuncs;

#define EXON_DATA \
  ANNOTATEDSEQFEATURE_DATA \
  int stickyRank; \
  FeatureSet supportingFeatures; \
  char *seqCacheString;

#define FUNCSTRUCTTYPE ExonFuncs
struct ExonStruct {
  EXON_DATA
};
#undef FUNCSTRUCTTYPE

#define Exon_setStart(exon,start) AnnotatedSeqFeature_setStart((exon),(start))
#define Exon_getStart(exon) AnnotatedSeqFeature_getStart((exon))

#define Exon_setEnd(exon,end) AnnotatedSeqFeature_setEnd((exon),(end))
#define Exon_getEnd(exon) AnnotatedSeqFeature_getEnd((exon))

#define Exon_setSeqCacheString(exon,seq) (exon)->seqCacheString = (seq)
#define Exon_getSeqCacheString(exon) (exon)->seqCacheString

#define Exon_setScore(exon,score) AnnotatedSeqFeature_setScore((exon),(score))
#define Exon_getScore(exon) AnnotatedSeqFeature_getScore((exon))

#define Exon_setpValue(exon,pValue) AnnotatedSeqFeature_setpValue((exon),(pValue))
#define Exon_getpValue(exon) AnnotatedSeqFeature_getpValue((exon))

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

#define Exon_getSupportingFeatureAt(exon,ind) FeatureSet_getFeatureAt(&((exon)->supportingFeatures),(ind))
#define Exon_addSupportingFeature(exon, sf) FeatureSet_addFeature(&((exon)->supportingFeatures),(sf))
#define Exon_getAllSupportingFeatures(exon) FeatureSet_getFeatures(&((exon)->supportingFeatures))
#define Exon_getSupportingFeatureCount(exon) FeatureSet_getNumFeature(&((exon)->supportingFeatures))


#define Exon_setContig(exon,c) AnnotatedSeqFeature_setContig((exon),(c))
#define Exon_getContig(exon) AnnotatedSeqFeature_getContig((exon))

Exon *Exon_transformRawContigToSlice(Exon *exon, Slice *slice);
Exon *Exon_transformSliceToRawContig(Exon *exon);

#define Exon_transformToSlice(exon,slice) SeqFeature_transformToSlice((exon), (slice))



Exon *Exon_new();
Exon *Exon_copy(Exon *copy, Exon *orig, CopyDepth depth);
char *Exon_getSeqString(Exon *exon);

int Exon_reverseStrandCompFunc(const void *a, const void *b);
int Exon_forwardStrandCompFunc(const void *a, const void *b);

void Exon_loadGenomicMapper(Exon *exon, Mapper *mapper, IDType id, int start);
Exon *Exon_adjustStartEnd(Exon *exon, int startAdjust, int endAdjust);


#ifdef __EXON_MAIN__
  ExonFuncs 
    exonFuncs = {
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
                 (Exon_TransformToRawContigFunc)SeqFeature_transformToRawContig,
                 (Exon_TransformToSliceFunc)SeqFeature_transformToSlice,
                 Exon_transformRawContigToSlice,
                 Exon_transformSliceToRawContig,
                 NULL // transformSliceToSlice
                };
#else
  extern ExonFuncs exonFuncs;
#endif

#endif
