#ifndef __BASEALIGNFEATURE_H__
#define __BASEALIGNFEATURE_H__

#include "DataModelTypes.h"
#include "SeqFeature.h"

struct BaseAlignFeatureStruct {
  SeqFeature sf1;
  SeqFeature sf2;
  char *cigarString;
  char *hitId;
};

BaseAlignFeature *BaseAlignFeature_new(void);

char *BaseAlignFeature_setCigarString(BaseAlignFeature *fp, char *ciggy);
#define BaseAlignFeature_getCigarString(fp)  (fp)->cigarString

#define BaseAlignFeature_setStableId(fp,stableId)  SeqFeature_setStableId(&((fp)->sf1),(stableId))
#define BaseAlignFeature_getStableId(fp)  SeqFeature_getStableId(&((fp)->sf1))

char *BaseAlignFeature_setCigarString(BaseAlignFeature *fp, char *hid);

#define BaseAlignFeature_getHitId(fp)  (fp)->hitId

#define BaseAlignFeature_setStart(fp,start) SeqFeature_setStart(&((fp)->sf1),(start))
#define BaseAlignFeature_getStart(fp) SeqFeature_getStart(&((fp)->sf1))

#define BaseAlignFeature_setHitStart(fp,start) SeqFeature_setStart(&((fp)->sf2),(start))
#define BaseAlignFeature_getHitStart(fp) SeqFeature_getStart(&((fp)->sf2))

#define BaseAlignFeature_setEnd(fp,end) SeqFeature_setEnd(&((fp)->sf1),(end))
#define BaseAlignFeature_getEnd(fp) SeqFeature_getEnd(&((fp)->sf1))

#define BaseAlignFeature_setHitEnd(fp,end) SeqFeature_setEnd(&((fp)->sf2),(end))
#define BaseAlignFeature_getHitEnd(fp) SeqFeature_getEnd(&((fp)->sf2))

#define BaseAlignFeature_setStrand(fp,strand) SeqFeature_setStrand(&((fp)->sf1),(strand))
#define BaseAlignFeature_getStrand(fp) SeqFeature_getStrand(&((fp)->sf1))

#define BaseAlignFeature_setHitStrand(fp,strand) SeqFeature_setStrand(&((fp)->sf2),(strand))
#define BaseAlignFeature_getHitStrand(fp) SeqFeature_getStrand(&((fp)->sf2))

#define BaseAlignFeature_setDbID(fp,dbID) SeqFeature_setDbID(&((fp)->sf1),(dbID))
#define BaseAlignFeature_getDbID(fp) SeqFeature_getDbID(&((fp)->sf1))

#define BaseAlignFeature_setAnalysis(fp,anal) SeqFeature_setAnalysis(&((fp)->sf1),(anal))
#define BaseAlignFeature_getAnalysis(fp) SeqFeature_getAnalysis(&((fp)->sf1))

#define BaseAlignFeature_setContig(fp,contig) SeqFeature_setContig(&((fp)->sf1),(contig))
#define BaseAlignFeature_getContig(fp) SeqFeature_getContig(&((fp)->sf1))

#define BaseAlignFeature_setScore(fp,score) SeqFeature_setScore(&((fp)->sf1),(score))
#define BaseAlignFeature_getScore(fp) SeqFeature_getScore(&((fp)->sf1))

#define BaseAlignFeature_setEValue(fp,ev) SeqFeature_setEValue(&((fp)->sf1),(ev))
#define BaseAlignFeature_getEValue(fp) SeqFeature_getEValue(&((fp)->sf1))

#define BaseAlignFeature_setPercId(fp,pid) SeqFeature_setPercId(&((fp)->sf1),(pid))
#define BaseAlignFeature_getPercId(fp) SeqFeature_getPercId(&((fp)->sf1))

#endif
