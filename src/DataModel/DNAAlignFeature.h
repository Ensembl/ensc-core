#ifndef __DNAALIGNFEATURE_H__
#define __DNAALIGNFEATURE_H__

#include "DataModelTypes.h"
#include "BaseAlignFeature.h"

struct DNAAlignFeatureStruct {
  BaseAlignFeature baf;
};

DNAAlignFeature *DNAAlignFeature_new(void);

#define DNAAlignFeature_setCigarString(fp, ciggy) BaseAlignFeature_setCigarString((BaseAlignFeature *)&((fp)->baf), (ciggy))
#define DNAAlignFeature_getCigarString(fp) BaseAlignFeature_getCigarString(&((fp)->baf))

#define DNAAlignFeature_setHitId(fp,stableId)  BaseAlignFeature_setHitId(&((fp)->baf),(stableId))
#define DNAAlignFeature_getHitId(fp)  BaseAlignFeature_getHitId(&((fp)->baf))

#define DNAAlignFeature_setStart(fp,start) BaseAlignFeature_setStart(&((fp)->baf),(start))
#define DNAAlignFeature_getStart(fp) BaseAlignFeature_getStart(&((fp)->baf))

#define DNAAlignFeature_setHitStart(fp,start) BaseAlignFeature_setHitStart(&((fp)->baf),(start))
#define DNAAlignFeature_getHitStart(fp) BaseAlignFeature_getHitStart(&((fp)->baf))

#define DNAAlignFeature_setEnd(fp,end) BaseAlignFeature_setEnd(&((fp)->baf),(end))
#define DNAAlignFeature_getEnd(fp) BaseAlignFeature_getEnd(&((fp)->baf))

#define DNAAlignFeature_setHitEnd(fp,end) BaseAlignFeature_setHitEnd(&((fp)->baf),(end))
#define DNAAlignFeature_getHitEnd(fp) BaseAlignFeature_getHitEnd(&((fp)->baf))

#define DNAAlignFeature_setStrand(fp,strand) BaseAlignFeature_setStrand(&((fp)->baf),(strand))
#define DNAAlignFeature_getStrand(fp) BaseAlignFeature_getStrand(&((fp)->baf))

#define DNAAlignFeature_setHitStrand(fp,strand) BaseAlignFeature_setStrand(&((fp)->baf),(strand))
#define DNAAlignFeature_getHitStrand(fp) BaseAlignFeature_getStrand(&((fp)->baf))

#define DNAAlignFeature_setDbID(fp,dbID) BaseAlignFeature_setDbID(&((fp)->baf),(dbID))
#define DNAAlignFeature_getDbID(fp) BaseAlignFeature_getDbID(&((fp)->baf))

#define DNAAlignFeature_setAnalysis(fp,anal) BaseAlignFeature_setAnalysis(&((fp)->baf),(anal))
#define DNAAlignFeature_getAnalysis(fp) BaseAlignFeature_getAnalysis(&((fp)->baf))

#define DNAAlignFeature_setContig(fp,contig) BaseAlignFeature_setContig(&((fp)->baf),(contig))
#define DNAAlignFeature_getContig(fp) BaseAlignFeature_getContig(&((fp)->baf))

#define DNAAlignFeature_setScore(fp,score) BaseAlignFeature_setScore(&((fp)->baf),(score))
#define DNAAlignFeature_getScore(fp) BaseAlignFeature_getScore(&((fp)->baf))

#define DNAAlignFeature_setEValue(fp,ev) BaseAlignFeature_setEValue(&((fp)->baf),(ev))
#define DNAAlignFeature_getEValue(fp) BaseAlignFeature_getEValue(&((fp)->baf))

#define DNAAlignFeature_setPercId(fp,pid) BaseAlignFeature_setPercId(&((fp)->baf),(pid))
#define DNAAlignFeature_getPercId(fp) BaseAlignFeature_getPercId(&((fp)->baf))


#endif
