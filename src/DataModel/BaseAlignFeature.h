#ifndef __BASEALIGNFEATURE_H__
#define __BASEALIGNFEATURE_H__

#include "DataModelTypes.h"
#include "FeaturePair.h"

#include "EnsRoot.h"

#define BASEALIGNFEATUREFUNC_TYPES(CLASSTYPE) \
 SEQFEATUREFUNC_TYPES(CLASSTYPE) \
 typedef int (* CLASSTYPE ## _GetHitUnitFunc)(void); \
 typedef int (* CLASSTYPE ## _GetQueryUnitFunc)(void);

#define BASEALIGNFEATUREFUNCS_DATA(CLASSTYPE) \
  FEATUREPAIRFUNCS_DATA(CLASSTYPE) \
  CLASSTYPE ## _GetHitUnitFunc getHitUnit; \
  CLASSTYPE ## _GetQueryUnitFunc getQueryUnit;


BASEALIGNFEATUREFUNC_TYPES(BaseAlignFeature)

typedef struct BaseAlignFeatureFuncsStruct {
  BASEALIGNFEATUREFUNCS_DATA(BaseAlignFeature)
} BaseAlignFeatureFuncs;
  

#define BASEALIGNFEATURE_DATA \
  FEATUREPAIR_DATA \
  char       *cigarString;


#define FUNCSTRUCTTYPE BaseAlignFeatureFuncs
struct BaseAlignFeatureStruct {
  BASEALIGNFEATURE_DATA
};
#undef FUNCSTRUCTTYPE

BaseAlignFeature *BaseAlignFeature_new(void);

char *BaseAlignFeature_setCigarString(BaseAlignFeature *fp, char *ciggy);
#define BaseAlignFeature_getCigarString(fp)  (fp)->cigarString

#define BaseAlignFeature_getHitSeqName(fp)  FeaturePair_getHitSeqName((fp))
#define BaseAlignFeature_setHitSeqName(fp,hid) FeaturePair_setHitSeqName((FeaturePair *)(fp),(hid))

#define BaseAlignFeature_setStart(fp,start) FeaturePair_setStart((fp),(start))
#define BaseAlignFeature_getStart(fp) FeaturePair_getStart((fp))

#define BaseAlignFeature_setHitStart(fp,start) FeaturePair_setHitStart((fp),(start))
#define BaseAlignFeature_getHitStart(fp) FeaturePair_getHitStart((fp))

#define BaseAlignFeature_setEnd(fp,end) FeaturePair_setEnd((fp),(end))
#define BaseAlignFeature_getEnd(fp) FeaturePair_getEnd((fp))

#define BaseAlignFeature_setHitEnd(fp,end) FeaturePair_setHitEnd((fp),(end))
#define BaseAlignFeature_getHitEnd(fp) FeaturePair_getHitEnd((fp))

#define BaseAlignFeature_setStrand(fp,strand) FeaturePair_setStrand((fp),(strand))
#define BaseAlignFeature_getStrand(fp) FeaturePair_getStrand((fp))

#define BaseAlignFeature_setEndPhase(fp,end) FeaturePair_setEndPhase((fp),(end))
#define BaseAlignFeature_getEndPhase(fp) FeaturePair_getEndPhase((fp))

#define BaseAlignFeature_setPhase(fp,end) FeaturePair_setPhase((fp),(end))
#define BaseAlignFeature_getPhase(fp) FeaturePair_getPhase((fp))

#define BaseAlignFeature_setHitStrand(fp,strand) FeaturePair_setHitStrand((fp),(strand))
#define BaseAlignFeature_getHitStrand(fp) FeaturePair_getHitStrand((fp))

#define BaseAlignFeature_setDbID(fp,dbID) FeaturePair_setDbID((fp),(dbID))
#define BaseAlignFeature_getDbID(fp) FeaturePair_getDbID((fp))

#define BaseAlignFeature_setAnalysis(fp,anal) FeaturePair_setAnalysis((fp),(anal))
#define BaseAlignFeature_getAnalysis(fp) FeaturePair_getAnalysis((fp))

#define BaseAlignFeature_setContig(fp,contig) FeaturePair_setContig((fp),(contig))
#define BaseAlignFeature_getContig(fp) FeaturePair_getContig((fp))

#define BaseAlignFeature_setScore(fp,score) FeaturePair_setScore((fp),(score))
#define BaseAlignFeature_getScore(fp) FeaturePair_getScore((fp))

#define BaseAlignFeature_setpValue(fp,ev) FeaturePair_setpValue((fp),(ev))
#define BaseAlignFeature_getpValue(fp) FeaturePair_getpValue((fp))

#define BaseAlignFeature_setPercId(fp,pid) FeaturePair_setPercId((fp),(pid))
#define BaseAlignFeature_getPercId(fp) FeaturePair_getPercId((fp))

#define BaseAlignFeature_setSeqName(fp,str) FeaturePair_setSeqName((fp),(str))

#define BaseAlignFeature_transformToRawContig(fp) FeaturePair_transformToRawContig((fp))
#define BaseAlignFeature_transformToSlice(fp,slice) FeaturePair_transformToSlice((fp),(slice))

Vector *BaseAlignFeature_parseCigar(BaseAlignFeature *baf);
Vector *BaseAlignFeature_transformSliceToRawContigImpl(BaseAlignFeature *baf);
int BaseAlignFeature_getHitUnitImpl(void);
int BaseAlignFeature_getQueryUnitImpl(void);
Vector *BaseAlignFeature_transformFeatureSliceToRawContig(BaseAlignFeature *baf, FeaturePair *fp);
Vector *BaseAlignFeature_getUngappedFeatures(BaseAlignFeature *baf);
void BaseAlignFeature_reverseComplementImpl(BaseAlignFeature *baf);
int BaseAlignFeature_parseFeatures(BaseAlignFeature *baf, Vector *features);

#define BaseAlignFeature_getHitUnit() \
      ((baf)->funcs->getHitUnit == NULL ? \
         (fprintf(stderr,"Error: Null pointer for getHitUnit - bye\n"),  exit(1), 0) : \
         ((baf)->funcs->getHitUnit()))

#define BaseAlignFeature_getQueryUnit() \
      ((baf)->funcs->getQueryUnit == NULL ? \
         (fprintf(stderr,"Error: Null pointer for getQueryUnit - bye\n"),  exit(1), 0) : \
         ((baf)->funcs->getQueryUnit()))



#ifdef __BASEALIGNFEATURE_MAIN__
  BaseAlignFeatureFuncs 
    baseAlignFeatureFuncs = {
                             NULL, // getStart
                             NULL, // setStart
                             NULL, // getEnd
                             NULL, // setEnd
                             NULL, // getStrand
                             NULL, // setStrand
                             NULL, // getSeq
                             NULL, // setSeq
                             NULL, // getLength
                             BaseAlignFeature_reverseComplementImpl,
                             (BaseAlignFeature_TransformToRawContigFunc)SeqFeature_transformToRawContigImpl,
                             (BaseAlignFeature_TransformToSliceFunc)SeqFeature_transformToSliceImpl,
                             (BaseAlignFeature_TransformRawContigToSliceFunc)SeqFeature_transformRawContigToSliceImpl, // Que???
                             BaseAlignFeature_transformSliceToRawContigImpl,
                             NULL, // transformSliceToSlice
                             BaseAlignFeature_getHitUnitImpl,
                             BaseAlignFeature_getQueryUnitImpl
                            };
#else 
  extern BaseAlignFeatureFuncs baseAlignFeatureFuncs;
#endif

#endif
