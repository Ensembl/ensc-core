#ifndef __FEATUREPAIR_H__
#define __FEATUREPAIR_H__

#include "DataModelTypes.h"
#include "SeqFeature.h"

#define FEATUREPAIRFUNCS_DATA(CLASSTYPE) \
  SEQFEATUREFUNCS_DATA(CLASSTYPE)

SEQFEATUREFUNC_TYPES(FeaturePair)

typedef struct FeaturePairFuncsStruct {
  FEATUREPAIRFUNCS_DATA(FeaturePair)
} FeaturePairFuncs;
  

#define FEATUREPAIR_DATA \
  SEQFEATURE_DATA \
  int         hitStart; \
  int         hitEnd; \
  signed char hitStrand; \
  signed char hitPhase; \
  signed char hitEndPhase; \
  char       *hitId;


#define FUNCSTRUCTTYPE FeaturePairFuncs
struct FeaturePairStruct {
  FEATUREPAIR_DATA
};
#undef FUNCSTRUCTTYPE

FeaturePair *FeaturePair_new(void);

#define FeaturePair_getHitSeqName(fp)  (fp)->hitId
char *FeaturePair_setHitSeqName(FeaturePair *fp, char *str);

#define FeaturePair_getSeqName(fp)  SeqFeature_getSeqName((fp))
#define FeaturePair_setSeqName(fp,str)  SeqFeature_setSeqName((fp),(str))

#define FeaturePair_setStart(fp,start) SeqFeature_setStart((fp),(start))
#define FeaturePair_getStart(fp) SeqFeature_getStart((fp))

#define FeaturePair_setHitStart(fp,start) (fp)->hitStart = start
#define FeaturePair_getHitStart(fp) (fp)->hitStart

#define FeaturePair_setEnd(fp,end) SeqFeature_setEnd((fp),(end))
#define FeaturePair_getEnd(fp) SeqFeature_getEnd((fp))

#define FeaturePair_setPhase(fp,phase) SeqFeature_setPhase((fp),(phase))
#define FeaturePair_getPhase(fp) SeqFeature_getPhase((fp))

#define FeaturePair_setEndPhase(fp,endPhase) SeqFeature_setEndPhase((fp),(endPhase))
#define FeaturePair_getEndPhase(fp) SeqFeature_getEndPhase((fp))

#define FeaturePair_setHitEnd(fp,end) (fp)->hitEnd = end
#define FeaturePair_getHitEnd(fp) (fp)->hitEnd

#define FeaturePair_setStrand(fp,strand) SeqFeature_setStrand((fp),(strand))
#define FeaturePair_getStrand(fp) SeqFeature_getStrand((fp))

#define FeaturePair_setHitStrand(fp,strand) (fp)->hitStrand = (strand)
#define FeaturePair_getHitStrand(fp) (fp)->hitStrand

#define FeaturePair_setDbID(fp,dbID) SeqFeature_setDbID((fp),(dbID))
#define FeaturePair_getDbID(fp) SeqFeature_getDbID((fp))

#define FeaturePair_setAnalysis(fp,anal) SeqFeature_setAnalysis((fp),(anal))
#define FeaturePair_getAnalysis(fp) SeqFeature_getAnalysis((fp))

#define FeaturePair_setContig(fp,contig) SeqFeature_setContig((fp),(contig))
#define FeaturePair_getContig(fp) SeqFeature_getContig((fp))

#define FeaturePair_setScore(fp,score) SeqFeature_setScore((fp),(score))
#define FeaturePair_getScore(fp) SeqFeature_getScore((fp))

#define FeaturePair_setpValue(fp,ev) SeqFeature_setpValue((fp),(ev))
#define FeaturePair_getpValue(fp) SeqFeature_getpValue((fp))

#define FeaturePair_setPercId(fp,pid) SeqFeature_setPercId((fp),(pid))
#define FeaturePair_getPercId(fp) SeqFeature_getPercId((fp))

#define FeaturePair_transformToRawContig(fp) SeqFeature_transformToRawContig((fp))
#define FeaturePair_transformToSlice(fp,slice) SeqFeature_transformToSlice((fp),(slice))

#ifdef __FEATUREPAIR_MAIN__
  FeaturePairFuncs 
    featurePairFuncs = {
                        NULL, // getStart
                        NULL, // setStart
                        NULL, // getEnd
                        NULL  // setEnd
                       };
#else 
  extern FeaturePairFuncs featurePairFuncs;
#endif

#endif
