#ifndef __DNAALIGNFEATURE_H__
#define __DNAALIGNFEATURE_H__

#include "DataModelTypes.h"
#include "BaseAlignFeature.h"

BASEALIGNFEATUREFUNC_TYPES(DNAAlignFeature)

typedef struct DNAAlignFeatureFuncsStruct {
  BASEALIGNFEATUREFUNCS_DATA(DNAAlignFeature)
} DNAAlignFeatureFuncs;
  


#define FUNCSTRUCTTYPE DNAAlignFeatureFuncs
struct DNAAlignFeatureStruct {
  BASEALIGNFEATURE_DATA
};
#undef FUNCSTRUCTTYPE

DNAAlignFeature *DNAAlignFeature_new(void);

#define DNAAlignFeature_setCigarString(fp, ciggy) BaseAlignFeature_setCigarString((BaseAlignFeature *)(fp), (ciggy))
#define DNAAlignFeature_getCigarString(fp) BaseAlignFeature_getCigarString((fp))

#define DNAAlignFeature_getUngappedFeatures(fp) BaseAlignFeature_getUngappedFeatures((fp))

#define DNAAlignFeature_setHitSeqName(fp,stableId)  BaseAlignFeature_setHitSeqName((BaseAlignFeature *)(fp),(stableId))
#define DNAAlignFeature_getHitSeqName(fp)  BaseAlignFeature_getHitSeqName((fp))

#define DNAAlignFeature_setStart(fp,start) BaseAlignFeature_setStart((fp),(start))
#define DNAAlignFeature_getStart(fp) BaseAlignFeature_getStart((fp))

#define DNAAlignFeature_setHitStart(fp,start) BaseAlignFeature_setHitStart((fp),(start))
#define DNAAlignFeature_getHitStart(fp) BaseAlignFeature_getHitStart((fp))

#define DNAAlignFeature_setEnd(fp,end) BaseAlignFeature_setEnd((fp),(end))
#define DNAAlignFeature_getEnd(fp) BaseAlignFeature_getEnd((fp))

#define DNAAlignFeature_setHitEnd(fp,end) BaseAlignFeature_setHitEnd((fp),(end))
#define DNAAlignFeature_getHitEnd(fp) BaseAlignFeature_getHitEnd((fp))

#define DNAAlignFeature_setStrand(fp,strand) BaseAlignFeature_setStrand((fp),(strand))
#define DNAAlignFeature_getStrand(fp) BaseAlignFeature_getStrand((fp))

#define DNAAlignFeature_setHitStrand(fp,strand) BaseAlignFeature_setHitStrand((fp),(strand))
#define DNAAlignFeature_getHitStrand(fp) BaseAlignFeature_getHitStrand((fp))

#define DNAAlignFeature_setDbID(fp,dbID) BaseAlignFeature_setDbID((fp),(dbID))
#define DNAAlignFeature_getDbID(fp) BaseAlignFeature_getDbID((fp))

#define DNAAlignFeature_setAnalysis(fp,anal) BaseAlignFeature_setAnalysis((fp),(anal))
#define DNAAlignFeature_getAnalysis(fp) BaseAlignFeature_getAnalysis((fp))

#define DNAAlignFeature_setContig(fp,contig) BaseAlignFeature_setContig((fp),(contig))
#define DNAAlignFeature_getContig(fp) BaseAlignFeature_getContig((fp))

#define DNAAlignFeature_setScore(fp,score) BaseAlignFeature_setScore((fp),(score))
#define DNAAlignFeature_getScore(fp) BaseAlignFeature_getScore((fp))

#define DNAAlignFeature_setpValue(fp,ev) BaseAlignFeature_setpValue((fp),(ev))
#define DNAAlignFeature_getpValue(fp) BaseAlignFeature_getpValue((fp))

#define DNAAlignFeature_setPercId(fp,pid) BaseAlignFeature_setPercId((fp),(pid))
#define DNAAlignFeature_getPercId(fp) BaseAlignFeature_getPercId((fp))

int DNAAlignFeature_getHitUnit(void);
int DNAAlignFeature_getQueryUnit(void);

#define DNAAlignFeature_transformToRawContig(fp) BaseAlignFeature_transformToRawContig((fp))
#define DNAAlignFeature_transformToSlice(fp,slice) BaseAlignFeature_transformToSlice((fp),(slice))

#ifdef __DNAALIGNFEATURE_MAIN__
  DNAAlignFeatureFuncs 
    dnaAlignFeatureFuncs = {
                             NULL, // getStart
                             NULL, // setStart
                             NULL, // getEnd
                             NULL, // setEnd
                             NULL, // getStrand
                             NULL, // setStrand
                             NULL, // getSeq
                             NULL, // setSeq
                             NULL, // getLength
                             (DNAAlignFeature_ReverseComplementFunc)BaseAlignFeature_reverseComplementImpl,
                             (DNAAlignFeature_TransformToRawContigFunc)SeqFeature_transformToRawContigImpl,
                             (DNAAlignFeature_TransformToSliceFunc)SeqFeature_transformToSliceImpl,
                             (DNAAlignFeature_TransformRawContigToSliceFunc)SeqFeature_transformRawContigToSliceImpl,
                             (DNAAlignFeature_TransformSliceToRawContigFunc)BaseAlignFeature_transformSliceToRawContigImpl,
                             NULL, // transformSliceToSlice
                             DNAAlignFeature_getHitUnit,
                             DNAAlignFeature_getQueryUnit
                            };
#else 
  extern DNAAlignFeatureFuncs dnaAlignFeatureFuncs;
#endif

#endif
