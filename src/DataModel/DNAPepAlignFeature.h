#ifndef __DNAPEPALIGNFEATURE_H__
#define __DNAPEPALIGNFEATURE_H__

#include "DataModelTypes.h"
#include "BaseAlignFeature.h"

struct DNAPepAlignFeatureStruct {
  BaseAlignFeature baf;
};

DNAPepAlignFeature *DNAPepAlignFeature_new(void);

#define DNAPepAlignFeature_setCigarString(fp, ciggy) BaseAlignFeature_setCigarString((BaseAlignFeature *)&((fp)->baf), (ciggy))
#define DNAPepAlignFeature_getCigarString(fp) BaseAlignFeature_getCigarString(&((fp)->baf))

#define DNAPepAlignFeature_setHitId(fp,stableId)  BaseAlignFeature_setHitId(&((fp)->baf),(stableId))
#define DNAPepAlignFeature_getHitId(fp)  BaseAlignFeature_getHitId(&((fp)->baf))

#define DNAPepAlignFeature_setStart(fp,start) BaseAlignFeature_setStart(&((fp)->baf),(start))
#define DNAPepAlignFeature_getStart(fp) BaseAlignFeature_getStart(&((fp)->baf))

#define DNAPepAlignFeature_setHitStart(fp,start) BaseAlignFeature_setHitStart(&((fp)->baf),(start))
#define DNAPepAlignFeature_getHitStart(fp) BaseAlignFeature_getHitStart(&((fp)->baf))

#define DNAPepAlignFeature_setEnd(fp,end) BaseAlignFeature_setEnd(&((fp)->baf),(end))
#define DNAPepAlignFeature_getEnd(fp) BaseAlignFeature_getEnd(&((fp)->baf))

#define DNAPepAlignFeature_setHitEnd(fp,end) BaseAlignFeature_setHitEnd(&((fp)->baf),(end))
#define DNAPepAlignFeature_getHitEnd(fp) BaseAlignFeature_getHitEnd(&((fp)->baf))

#define DNAPepAlignFeature_setStrand(fp,strand) BaseAlignFeature_setStrand(&((fp)->baf),(strand))
#define DNAPepAlignFeature_getStrand(fp) BaseAlignFeature_getStrand(&((fp)->baf))

#define DNAPepAlignFeature_setHitStrand(fp,strand) BaseAlignFeature_setStrand(&((fp)->baf),(strand))
#define DNAPepAlignFeature_getHitStrand(fp) BaseAlignFeature_getStrand(&((fp)->baf))

#define DNAPepAlignFeature_setDbID(fp,dbID) BaseAlignFeature_setDbID(&((fp)->baf),(dbID))
#define DNAPepAlignFeature_getDbID(fp) BaseAlignFeature_getDbID(&((fp)->baf))

#define DNAPepAlignFeature_setAnalysis(fp,anal) BaseAlignFeature_setAnalysis(&((fp)->baf),(anal))
#define DNAPepAlignFeature_getAnalysis(fp) BaseAlignFeature_getAnalysis(&((fp)->baf))

#define DNAPepAlignFeature_setContig(fp,contig) BaseAlignFeature_setContig(&((fp)->baf),(contig))
#define DNAPepAlignFeature_getContig(fp) BaseAlignFeature_getContig(&((fp)->baf))

#define DNAPepAlignFeature_setScore(fp,score) BaseAlignFeature_setScore(&((fp)->baf),(score))
#define DNAPepAlignFeature_getScore(fp) BaseAlignFeature_getScore(&((fp)->baf))

#define DNAPepAlignFeature_setEValue(fp,ev) BaseAlignFeature_setEValue(&((fp)->baf),(ev))
#define DNAPepAlignFeature_getEValue(fp) BaseAlignFeature_getEValue(&((fp)->baf))

#define DNAPepAlignFeature_setPercId(fp,pid) BaseAlignFeature_setPercId(&((fp)->baf),(pid))
#define DNAPepAlignFeature_getPercId(fp) BaseAlignFeature_getPercId(&((fp)->baf))

#endif
