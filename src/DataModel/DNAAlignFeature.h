#ifndef __DNAALIGNFEATURE_H__
#define __DNAALIGNFEATURE_H__

#include "DataModelTypes.h"
#include "BaseAlignFeature.h"

struct DNAAlignFeatureStruct {
  SeqFeature sf1;
  SeqFeature sf2;
  char *cigarString;
};

DNAAlignFeature *DNAAlignFeature_new(void);

#define DNAAlignFeature_setCigarString(fp, ciggy) BaseAlignFeature_setCigarString((fp), (ciggy))
#define DNAAlignFeature_getCigarString(fp) BaseAlignFeature_getCigarString((fp));

#define DNAAlignFeature_setStableId(fp,stableId)  BaseAlignFeature_setStableId((fp),(stableId))
#define DNAAlignFeature_getStableId(fp)  BaseAlignFeature_getStableId((fp))

#define DNAAlignFeature_setHitId(fp,stableId)  BaseAlignFeature_setStableId((fp),(stableId))
#define DNAAlignFeature_getHitId(fp)  BaseAlignFeature_getStableId((fp))

#define DNAAlignFeature_setStart(fp,start) BaseAlignFeature_setStart((fp),(start))
#define DNAAlignFeature_getStart(fp) BaseAlignFeature_getStart((fp))

#define DNAAlignFeature_setHitStart(fp,start) BaseAlignFeature_setStart((fp),(start))
#define DNAAlignFeature_getHitStart(fp) BaseAlignFeature_getStart((fp))

#define DNAAlignFeature_setEnd(fp,end) BaseAlignFeature_setEnd((fp),(end))
#define DNAAlignFeature_getEnd(fp) BaseAlignFeature_getEnd((fp))

#define DNAAlignFeature_setHitEnd(fp,end) BaseAlignFeature_setEnd((fp),(end))
#define DNAAlignFeature_getHitEnd(fp) BaseAlignFeature_getEnd((fp))

#define DNAAlignFeature_setStrand(fp,strand) BaseAlignFeature_setStrand((fp),(strand))
#define DNAAlignFeature_getStrand(fp) BaseAlignFeature_getStrand((fp))

#define DNAAlignFeature_setHitStrand(fp,strand) BaseAlignFeature_setStrand((fp),(strand))
#define DNAAlignFeature_getHitStrand(fp) BaseAlignFeature_getStrand((fp))

#define DNAAlignFeature_setDbID(fp,dbID) BaseAlignFeature_setDbID((fp),(dbID))
#define DNAAlignFeature_getDbID(fp,dbID) BaseAlignFeature_getDbID((fp))

#endif
