#ifndef __BASECONTIG_H__
#define __BASECONTIG_H__

#include "DataModelTypes.h"
#include "Storable.h"

typedef enum ContigTypeEnum {
  CONTIGTYPE_NONE,
  RAWCONTIG,
  SLICE
} ContigType;

#define BASECONTIG_DATA \
  ContigType contigType; \
  Storable st; \
  int start; \
  int end;

struct BaseContigStruct {
  BASECONTIG_DATA
};

#define BaseContig_setDbID(bc,dbID) Storable_setDbID(&((bc)->st),dbID)
#define BaseContig_getDbID(bc) Storable_getDbID(&((bc)->st))

#define BaseContig_setAdaptor(bc,ad) Storable_setAdaptor(&((bc)->st),ad)
#define BaseContig_getAdaptor(bc) Storable_getAdaptor(&((bc)->st))

#define BaseContig_setContigType(bc,ct) (bc)->contigType = (ct)
#define BaseContig_getContigType(bc) (bc)->contigType

char *BaseContig_getName(BaseContig *bc);

#endif
