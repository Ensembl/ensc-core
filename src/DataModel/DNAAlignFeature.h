#ifndef __DNAALIGNFEATURE_H__
#define __DNAALIGNFEATURE_H__

#include "DataModelTypes.h"
#include "BaseAlignFeature.h"

#define FUNCSTRUCTTYPE BaseAlignFeatureFuncs
struct DNAAlignFeatureStruct {
  BASEALIGNFEATURE_DATA
};
#undef FUNCSTRUCTTYPE

DNAAlignFeature *DNAAlignFeature_new(void);

#define DNAAlignFeature_setCigarString(fp, ciggy) BaseAlignFeature_setCigarString((BaseAlignFeature *)(fp), (ciggy))
#define DNAAlignFeature_getCigarString(fp) BaseAlignFeature_getCigarString((fp))

#define DNAAlignFeature_setHitId(fp,stableId)  BaseAlignFeature_setHitId((BaseAlignFeature *)(fp),(stableId))
#define DNAAlignFeature_getHitId(fp)  BaseAlignFeature_getHitId((fp))

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

#define DNAAlignFeature_setHitStrand(fp,strand) BaseAlignFeature_setStrand((fp),(strand))
#define DNAAlignFeature_getHitStrand(fp) BaseAlignFeature_getStrand((fp))

#define DNAAlignFeature_setDbID(fp,dbID) BaseAlignFeature_setDbID((fp),(dbID))
#define DNAAlignFeature_getDbID(fp) BaseAlignFeature_getDbID((fp))

#define DNAAlignFeature_setAnalysis(fp,anal) BaseAlignFeature_setAnalysis((fp),(anal))
#define DNAAlignFeature_getAnalysis(fp) BaseAlignFeature_getAnalysis((fp))

#define DNAAlignFeature_setContig(fp,contig) BaseAlignFeature_setContig((fp),(contig))
#define DNAAlignFeature_getContig(fp) BaseAlignFeature_getContig((fp))

#define DNAAlignFeature_setScore(fp,score) BaseAlignFeature_setScore((fp),(score))
#define DNAAlignFeature_getScore(fp) BaseAlignFeature_getScore((fp))

#define DNAAlignFeature_setEValue(fp,ev) BaseAlignFeature_setEValue((fp),(ev))
#define DNAAlignFeature_getEValue(fp) BaseAlignFeature_getEValue((fp))

#define DNAAlignFeature_setPercId(fp,pid) BaseAlignFeature_setPercId((fp),(pid))
#define DNAAlignFeature_getPercId(fp) BaseAlignFeature_getPercId((fp))


#endif
