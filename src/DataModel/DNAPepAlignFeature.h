#ifndef __DNAPEPALIGNFEATURE_H__
#define __DNAPEPALIGNFEATURE_H__

#include "DataModelTypes.h"
#include "BaseAlignFeature.h"

struct DNAPepAlignFeatureStruct {
  SeqFeature sf1;
  SeqFeature sf2;
  char *cigarString;
};

DNAPepAlignFeature *DNAPepAlignFeature_new(void);

#define DNAPepAlignFeature_setCigarString(fp, ciggy) BaseAlignFeature_setCigarString((BaseAlignFeature *)(fp), (ciggy))
#define DNAPepAlignFeature_getCigarString(fp) BaseAlignFeature_getCigarString((fp));

#define DNAPepAlignFeature_setHitId(fp,stableId)  BaseAlignFeature_setHitId((fp),(stableId))
#define DNAPepAlignFeature_getHitId(fp)  BaseAlignFeature_getHitId((fp))

#define DNAPepAlignFeature_setStart(fp,start) BaseAlignFeature_setStart((fp),(start))
#define DNAPepAlignFeature_getStart(fp) BaseAlignFeature_getStart((fp))

#define DNAPepAlignFeature_setHitStart(fp,start) BaseAlignFeature_setStart((fp),(start))
#define DNAPepAlignFeature_getHitStart(fp) BaseAlignFeature_getStart((fp))

#define DNAPepAlignFeature_setEnd(fp,end) BaseAlignFeature_setEnd((fp),(end))
#define DNAPepAlignFeature_getEnd(fp) BaseAlignFeature_getEnd((fp))

#define DNAPepAlignFeature_setHitEnd(fp,end) BaseAlignFeature_setEnd((fp),(end))
#define DNAPepAlignFeature_getHitEnd(fp) BaseAlignFeature_getEnd((fp))

#define DNAPepAlignFeature_setStrand(fp,strand) BaseAlignFeature_setStrand((fp),(strand))
#define DNAPepAlignFeature_getStrand(fp) BaseAlignFeature_getStrand((fp))

#define DNAPepAlignFeature_setHitStrand(fp,strand) BaseAlignFeature_setStrand((fp),(strand))
#define DNAPepAlignFeature_getHitStrand(fp) BaseAlignFeature_getStrand((fp))

#define DNAPepAlignFeature_setDbID(fp,dbID) BaseAlignFeature_setDbID((fp),(dbID))
#define DNAPepAlignFeature_getDbID(fp) BaseAlignFeature_getDbID((fp))

#define DNAPepAlignFeature_setAnalysis(fp,anal) BaseAlignFeature_setAnalysis((fp),(anal))
#define DNAPepAlignFeature_getAnalysis(fp) BaseAlignFeature_getAnalysis((fp))

#define DNAPepAlignFeature_setContig(fp,contig) BaseAlignFeature_setContig((fp),(contig))
#define DNAPepAlignFeature_getContig(fp) BaseAlignFeature_getContig((fp))

#define DNAPepAlignFeature_setScore(fp,score) BaseAlignFeature_setScore((fp),(score))
#define DNAPepAlignFeature_getScore(fp) BaseAlignFeature_getScore((fp))

#define DNAPepAlignFeature_setEValue(fp,ev) BaseAlignFeature_setEValue((fp),(ev))
#define DNAPepAlignFeature_getEValue(fp) BaseAlignFeature_getEValue((fp))

#define DNAPepAlignFeature_setPercId(fp,pid) BaseAlignFeature_setPercId((fp),(pid))
#define DNAPepAlignFeature_getPercId(fp) BaseAlignFeature_getPercId((fp))

#endif
