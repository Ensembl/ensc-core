#ifndef __STABLEIDINFO_H__
#define __STABLEIDINFO_H__

#include <time.h>
#include <stdio.h>

typedef struct StableIdInfoStruct {
  char *stableId;
  int   version;
  time_t created;
  time_t modified;
} StableIdInfo;

char *StableIdInfo_setStableId(StableIdInfo *si, char *sid); 
#define StableIdInfo_getStableId(si)  (si)->stableId

#define StableIdInfo_setCreated(si,cd)  (si)->created = (cd)
#define StableIdInfo_getCreated(si)  (si)->created

#define StableIdInfo_setModified(si,mod)  (si)->modified = (mod)
#define StableIdInfo_getModified(si) (si)->modified

#define StableIdInfo_setVersion(si,ver)  (si)->version = (ver)
#define StableIdInfo_getVersion(si) (si)->version

#endif
