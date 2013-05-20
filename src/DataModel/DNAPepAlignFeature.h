#ifndef __DNAPEPALIGNFEATURE_H__
#define __DNAPEPALIGNFEATURE_H__

#include "DataModelTypes.h"
#include "BaseAlignFeature.h"

BASEALIGNFEATUREFUNC_TYPES(DNAPepAlignFeature)

typedef struct DNAPepAlignFeatureFuncsStruct {
  BASEALIGNFEATUREFUNCS_DATA(DNAPepAlignFeature)
} DNAPepAlignFeatureFuncs;
  
#define FUNCSTRUCTTYPE DNAPepAlignFeatureFuncs
struct DNAPepAlignFeatureStruct {
  BASEALIGNFEATURE_DATA
};
#undef FUNCSTRUCTTYPE

DNAPepAlignFeature *DNAPepAlignFeature_new(void);

#define DNAPepAlignFeature_setCigarString(fp, ciggy) BaseAlignFeature_setCigarString((BaseAlignFeature *)(fp), (ciggy))
#define DNAPepAlignFeature_getCigarString(fp) BaseAlignFeature_getCigarString((fp))

#define DNAPepAlignFeature_getUngappedFeatures(fp) BaseAlignFeature_getUngappedFeatures((fp))


#define DNAPepAlignFeature_setHitSeqName(fp,hname)  BaseAlignFeature_setHitSeqName((BaseAlignFeature *)(fp),(hname))
#define DNAPepAlignFeature_getHitSeqName(fp)  BaseAlignFeature_getHitSeqName((fp))

#define DNAPepAlignFeature_setStart(fp,start) BaseAlignFeature_setStart((fp),(start))
#define DNAPepAlignFeature_getStart(fp) BaseAlignFeature_getStart((fp))

#define DNAPepAlignFeature_setHitStart(fp,start) BaseAlignFeature_setHitStart((fp),(start))
#define DNAPepAlignFeature_getHitStart(fp) BaseAlignFeature_getHitStart((fp))

#define DNAPepAlignFeature_setEnd(fp,end) BaseAlignFeature_setEnd((fp),(end))
#define DNAPepAlignFeature_getEnd(fp) BaseAlignFeature_getEnd((fp))

#define DNAPepAlignFeature_setHitEnd(fp,end) BaseAlignFeature_setHitEnd((fp),(end))
#define DNAPepAlignFeature_getHitEnd(fp) BaseAlignFeature_getHitEnd((fp))

#define DNAPepAlignFeature_setStrand(fp,strand) BaseAlignFeature_setStrand((fp),(strand))
#define DNAPepAlignFeature_getStrand(fp) BaseAlignFeature_getStrand((fp))

#define DNAPepAlignFeature_setHitStrand(fp,strand) BaseAlignFeature_setHitStrand((fp),(strand))
#define DNAPepAlignFeature_getHitStrand(fp) BaseAlignFeature_getHitStrand((fp))

#define DNAPepAlignFeature_setDbID(fp,dbID) BaseAlignFeature_setDbID((fp),(dbID))
#define DNAPepAlignFeature_getDbID(fp) BaseAlignFeature_getDbID((fp))

#define DNAPepAlignFeature_setAnalysis(fp,anal) BaseAlignFeature_setAnalysis((fp),(anal))
#define DNAPepAlignFeature_getAnalysis(fp) BaseAlignFeature_getAnalysis((fp))

#define DNAPepAlignFeature_setContig(fp,contig) BaseAlignFeature_setContig((fp),(contig))
#define DNAPepAlignFeature_getContig(fp) BaseAlignFeature_getContig((fp))

#define DNAPepAlignFeature_setScore(fp,score) BaseAlignFeature_setScore((fp),(score))
#define DNAPepAlignFeature_getScore(fp) BaseAlignFeature_getScore((fp))

#define DNAPepAlignFeature_setpValue(fp,ev) BaseAlignFeature_setpValue((fp),(ev))
#define DNAPepAlignFeature_getpValue(fp) BaseAlignFeature_getpValue((fp))

#define DNAPepAlignFeature_setPercId(fp,pid) BaseAlignFeature_setPercId((fp),(pid))
#define DNAPepAlignFeature_getPercId(fp) BaseAlignFeature_getPercId((fp))

int DNAPepAlignFeature_getHitUnit(void);
int DNAPepAlignFeature_getQueryUnit(void);

#define DNAPepAlignFeature_transformToSlice(fp,slice) BaseAlignFeature_transformToSlice((fp),(slice))
#define DNAPepAlignFeature_transformToRawContig(fp) BaseAlignFeature_transformToRawContig((fp))

#define DNAPepAlignFeature_free(fp) BaseAlignFeature_free((fp))
#define DNAPepAlignFeature_shallowCopy(fp) DNAPepAlignFeature_shallowCopyImpl((fp))

void DNAPepAlignFeature_freeImpl(DNAPepAlignFeature *dpaf);
DNAPepAlignFeature *DNAPepAlignFeature_shallowCopyImpl(DNAPepAlignFeature *dpaf);


#ifdef __DNAPEPALIGNFEATURE_MAIN__
  DNAPepAlignFeatureFuncs
    dnaPepAlignFeatureFuncs = {
                             DNAPepAlignFeature_freeImpl,
                             DNAPepAlignFeature_shallowCopyImpl, // shallowCopy
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
                             (DNAPepAlignFeature_ReverseComplementFunc)BaseAlignFeature_reverseComplementImpl,
                             (DNAPepAlignFeature_TransformToRawContigFunc)SeqFeature_transformToRawContigImpl,
                             (DNAPepAlignFeature_TransformToSliceFunc)SeqFeature_transformToSliceImpl,
                             (DNAPepAlignFeature_TransformRawContigToSliceFunc)SeqFeature_transformRawContigToSliceImpl,
                             (DNAPepAlignFeature_TransformSliceToRawContigFunc)BaseAlignFeature_transformSliceToRawContigImpl,
                             NULL, // transformSliceToSlice
                             DNAPepAlignFeature_getHitUnit,
                             DNAPepAlignFeature_getQueryUnit
                            };
#else
  extern DNAPepAlignFeatureFuncs dnaPepAlignFeatureFuncs;
#endif

#endif
