#ifndef __EXON_H__
#define __EXON_H__

#include <time.h>

#include "EnsC.h"
#include "DataModelTypes.h"
#include "AnnotatedSeqFeature.h"
#include "StableIdInfo.h"
#include "Slice.h"
#include "SeqFeature.h"

#define EXONFUNC_TYPES(CLASSTYPE) \
  ANNOTATEDSEQFEATUREFUNC_TYPES(CLASSTYPE) \
  typedef void     (* CLASSTYPE ## _LoadGenomicMapperFunc)(CLASSTYPE *exon, Mapper *mapper, IDType id, int start); \
  typedef CLASSTYPE *   (* CLASSTYPE ## _AdjustStartEndFunc)(CLASSTYPE *exon, int startAdjust, int endAdjust); \
  typedef char *   (* CLASSTYPE ## _GetPeptideFunc)(CLASSTYPE *exon, Transcript *trans); \
  typedef void     (* CLASSTYPE ## _AddSupportingFeaturesFunc)(CLASSTYPE *exon, Vector *features); \
  typedef Vector * (* CLASSTYPE ## _GetAllSupportingFeaturesFunc)(CLASSTYPE *exon); \
  typedef char *   (* CLASSTYPE ## _GetSeqStringFunc)(CLASSTYPE *exon); \

EXONFUNC_TYPES(Exon)

#define EXONFUNCS_DATA(CLASSTYPE) \
  ANNOTATEDSEQFEATUREFUNCS_DATA(CLASSTYPE) \
  CLASSTYPE ## _LoadGenomicMapperFunc loadGenomicMapper; \
  CLASSTYPE ## _AdjustStartEndFunc adjustStartEnd; \
  CLASSTYPE ## _GetPeptideFunc getPeptide; \
  CLASSTYPE ## _AddSupportingFeaturesFunc addSupportingFeatures; \
  CLASSTYPE ## _GetAllSupportingFeaturesFunc getAllSupportingFeatures; \
  CLASSTYPE ## _GetSeqStringFunc getSeqString;

typedef struct ExonFuncsStruct {
  EXONFUNCS_DATA(Exon)
} ExonFuncs;

#define EXON_DATA \
  ANNOTATEDSEQFEATURE_DATA \
  char *seqCacheString; \
  Vector *supportingFeatures; \
  char isConstitutive; \

#define FUNCSTRUCTTYPE ExonFuncs
struct ExonStruct {
  EXON_DATA
};
#undef FUNCSTRUCTTYPE

#define Exon_isStored(exon, db) Storable_isStored(&((exon)->st), (db))

#define Exon_addFlag(exon, f) SeqFeature_addFlag((exon), (f))
#define Exon_getFlags(exon) SeqFeature_getFlags(exon)
#define Exon_removeFlag(exon, f) SeqFeature_removeFlag((exon), (f))

#define Exon_getSeqRegionStart(exon) SeqFeature_getSeqRegionStart((exon))
#define Exon_getSeqRegionEnd(exon) SeqFeature_getSeqRegionEnd((exon))
#define Exon_getSeqRegionStrand(exon) SeqFeature_getSeqRegionStrand((exon))

#define Exon_setIsCurrent(exon,isC)  StableIdInfo_setIsCurrent(&((exon)->si),(isC))
#define Exon_getIsCurrent(exon)  StableIdInfo_getIsCurrent(&((exon)->si))

#define Exon_setIsConstitutive(exon,isC) (exon)->isConstitutive = (isC)
#define Exon_getIsConstitutive(exon) (exon)->isConstitutive

#define Exon_setStart(exon,start) AnnotatedSeqFeature_setStart((exon),(start))
#define Exon_getStart(exon) AnnotatedSeqFeature_getStart((exon))

#define Exon_setEnd(exon,end) AnnotatedSeqFeature_setEnd((exon),(end))
#define Exon_getEnd(exon) AnnotatedSeqFeature_getEnd((exon))

//#define Exon_setSeqCacheString(exon,seq) (exon)->seqCacheString = (seq)
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

Vector *Exon_getAllSupportingFeaturesImpl(Exon *exon);
void Exon_addSupportingFeaturesImpl(Exon *exon, Vector *features);
void Exon_addSupportingFeature(Exon *exon, SeqFeature *sf);

void Exon_setSeqCacheString(Exon *exon, char *seq);


#define Exon_setContig(exon,c) AnnotatedSeqFeature_setContig((exon),(c))
#define Exon_getContig(exon) AnnotatedSeqFeature_getContig((exon))

#define Exon_setSlice(exon,sl) AnnotatedSeqFeature_setSlice((exon),(sl))
#define Exon_getSlice(exon) AnnotatedSeqFeature_getSlice((exon))


Exon *Exon_transformRawContigToSliceImpl(Exon *exon, Slice *slice);
Exon *Exon_transformSliceToRawContigImpl(Exon *exon);

#define Exon_transformToSlice(exon,slice) AnnotatedSeqFeature_transformToSlice((exon), (slice))
#define Exon_transformToRawContig(exon) AnnotatedSeqFeature_transformToRawContig((exon))

void Exon_flushSupportingFeatures(Exon *exon);

void Exon_freeImpl(Exon *exon);
#define Exon_free(exon) AnnotatedSeqFeature_free((exon))

Exon *Exon_shallowCopyImpl(Exon *exon);
#define Exon_shallowCopy(exon) AnnotatedSeqFeature_shallowCopy((exon))

Exon *Exon_new();
Exon *Exon_copy(Exon *copy, Exon *orig, CopyDepth depth);
char *Exon_getSeqStringImpl(Exon *exon);
char *Exon_getPeptideImpl(Exon *exon, Transcript *trans);

void Exon_getHashKey(Exon *exon, char *hashKey);

int Exon_reverseStrandCompFunc(const void *a, const void *b);
int Exon_forwardStrandCompFunc(const void *a, const void *b);

void Exon_loadGenomicMapperImpl(Exon *exon, Mapper *mapper, IDType id, int start);
Exon *Exon_adjustStartEndImpl(Exon *exon, int startAdjust, int endAdjust);

Exon *Exon_transfer(Exon *exon, Slice *slice);

#define Exon_loadGenomicMapper(exon,mapper,id,start) \
      ((exon)->funcs->loadGenomicMapper == NULL ? \
         (fprintf(stderr,"Error: Null pointer for loadGenomicMapper - bye\n"),  exit(1)) : \
         ((exon)->funcs->loadGenomicMapper((exon),(mapper),(id),(start))))

#define Exon_adjustStartEnd(exon,start,end) \
      ((exon)->funcs->adjustStartEnd == NULL ? \
         (fprintf(stderr,"Error: Null pointer for adjustStartEnd - bye\n"),  exit(1), (Exon *)NULL) : \
         ((exon)->funcs->adjustStartEnd((exon),(start),(end))))

#define Exon_addSupportingFeatures(exon,features) \
      ((exon)->funcs->addSupportingFeatures == NULL ? \
         (fprintf(stderr,"Error: Null pointer for addSupportingFeatures - bye\n"),  exit(1)) : \
         ((exon)->funcs->addSupportingFeatures((exon),(features))))

#define Exon_getAllSupportingFeatures(exon) \
      ((exon)->funcs->addSupportingFeatures == NULL ? \
         (fprintf(stderr,"Error: Null pointer for getAllSupportingFeatures - bye\n"),  exit(1), (Vector *)NULL) : \
         ((exon)->funcs->getAllSupportingFeatures((exon))))

#define Exon_getSeqString(exon) \
      ((exon)->funcs->getSeqString == NULL ? \
         (fprintf(stderr,"Error: Null pointer for getSeqString - bye\n"),  exit(1), (char *)NULL) : \
         ((exon)->funcs->getSeqString((exon))))

#define Exon_getPeptide(exon) \
      ((exon)->funcs->getPeptide == NULL ? \
         (fprintf(stderr,"Error: Null pointer for getPeptide - bye\n"),  exit(1), (char *)NULL) : \
         ((exon)->funcs->getPeptide((exon))))


#ifdef __EXON_MAIN__
  ExonFuncs 
    exonFuncs = {
                 Exon_freeImpl, // free
                 Exon_shallowCopyImpl, // shallowCopy
                 NULL, // deepCopy
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
                 (Exon_TransformToRawContigFunc)SeqFeature_transformToRawContigImpl,
                 (Exon_TransformToSliceFunc)SeqFeature_transformToSliceImpl,
                 Exon_transformRawContigToSliceImpl,
                 Exon_transformSliceToRawContigImpl,
                 NULL, // transformSliceToSlice
                 Exon_loadGenomicMapperImpl,
                 Exon_adjustStartEndImpl,
                 Exon_getPeptideImpl,
                 Exon_addSupportingFeaturesImpl,
                 Exon_getAllSupportingFeaturesImpl,
                 Exon_getSeqStringImpl
                };
#else
  extern ExonFuncs exonFuncs;
#endif

#endif
