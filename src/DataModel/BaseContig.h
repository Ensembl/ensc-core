#ifndef __BASECONTIG_H__
#define __BASECONTIG_H__

#include "DataModelTypes.h"
#include "Storable.h"
#include "EnsRoot.h"

typedef char * (*BaseContig_GetNameFunc)(BaseContig *bc);

#define BASECONTIGFUNCS_DATA \
  BaseContig_GetNameFunc *getName;
  
  
typedef struct BaseContigFuncsStruct {
  BASECONTIGFUNCS_DATA
} BaseContigFuncs;

#define BASECONTIG_DATA \
  ENSROOT_DATA \
  Storable st; \
  int start; \
  int end;

#undef FUNCSTRUCTTYPE
#define FUNCSTRUCTTYPE BaseContigFuncs
struct BaseContigStruct {
  BASECONTIG_DATA
};
#undef FUNCSTRUCTTYPE

#define BaseContig_setDbID(bc,dbID) Storable_setDbID(&((bc)->st),dbID)
#define BaseContig_getDbID(bc) Storable_getDbID(&((bc)->st))

#define BaseContig_setAdaptor(bc,ad) Storable_setAdaptor(&((bc)->st),ad)
#define BaseContig_getAdaptor(bc) Storable_getAdaptor(&((bc)->st))

char *BaseContig_getName(BaseContig *bc);

#define BaseContig_getObjectType(bc) (bc)->objectType

#endif
