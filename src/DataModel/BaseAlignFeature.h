#ifndef __BASEALIGNFEATURE_H__
#define __BASEALIGNFEATURE_H__

#include "DataModelTypes.h"
#include "SeqFeature.h"

#include "EnsRoot.h"

#define BASEALIGNFEATUREFUNCS_DATA \
  SEQFEATUREFUNCS_DATA \

typedef struct BaseAlignFeatureFuncsStruct {
  BASEALIGNFEATUREFUNCS_DATA
} BaseAlignFeatureFuncs;
  

#define BASEALIGNFEATURE_DATA \
  SEQFEATURE_DATA \
  int         hitStart; \
  int         hitEnd; \
  signed char hitPhase; \
  signed char hitEndPhase; \
  signed char hitStrand; \
  char       *cigarString; \
  char       *hitId;


#define FUNCSTRUCTTYPE BaseAlignFeatureFuncs
struct BaseAlignFeatureStruct {
  BASEALIGNFEATURE_DATA
};
#undef FUNCSTRUCTTYPE

BaseAlignFeature *BaseAlignFeature_new(void);

char *BaseAlignFeature_setCigarString(BaseAlignFeature *fp, char *ciggy);
#define BaseAlignFeature_getCigarString(fp)  (fp)->cigarString

char *BaseAlignFeature_setCigarString(BaseAlignFeature *fp, char *hid);

#define BaseAlignFeature_getHitId(fp)  (fp)->hitId
char *BaseAlignFeature_setHitId(BaseAlignFeature *baf, char *str);

#define BaseAlignFeature_setStart(fp,start) SeqFeature_setStart((fp),(start))
#define BaseAlignFeature_getStart(fp) SeqFeature_getStart((fp))

#define BaseAlignFeature_setHitStart(fp,start) (fp)->hitStart = start
#define BaseAlignFeature_getHitStart(fp) (fp)->hitStart

#define BaseAlignFeature_setEnd(fp,end) SeqFeature_setEnd((fp),(end))
#define BaseAlignFeature_getEnd(fp) SeqFeature_getEnd((fp))

#define BaseAlignFeature_setHitEnd(fp,end) (fp)->hitEnd = end
#define BaseAlignFeature_getHitEnd(fp) (fp)->hitEnd

#define BaseAlignFeature_setStrand(fp,strand) SeqFeature_setStrand((fp),(strand))
#define BaseAlignFeature_getStrand(fp) SeqFeature_getStrand((fp))

#define BaseAlignFeature_setHitStrand(fp,strand) (fp)->hitStrand = strand
#define BaseAlignFeature_getHitStrand(fp) (fp)->hitStrand

#define BaseAlignFeature_setDbID(fp,dbID) SeqFeature_setDbID((fp),(dbID))
#define BaseAlignFeature_getDbID(fp) SeqFeature_getDbID((fp))

#define BaseAlignFeature_setAnalysis(fp,anal) SeqFeature_setAnalysis((fp),(anal))
#define BaseAlignFeature_getAnalysis(fp) SeqFeature_getAnalysis((fp))

#define BaseAlignFeature_setContig(fp,contig) SeqFeature_setContig((fp),(contig))
#define BaseAlignFeature_getContig(fp) SeqFeature_getContig((fp))

#define BaseAlignFeature_setScore(fp,score) SeqFeature_setScore((fp),(score))
#define BaseAlignFeature_getScore(fp) SeqFeature_getScore((fp))

#define BaseAlignFeature_setEValue(fp,ev) SeqFeature_setEValue((fp),(ev))
#define BaseAlignFeature_getEValue(fp) SeqFeature_getEValue((fp))

#define BaseAlignFeature_setPercId(fp,pid) SeqFeature_setPercId((fp),(pid))
#define BaseAlignFeature_getPercId(fp) SeqFeature_getPercId((fp))

#endif
