#ifndef __SEQREGIONENTRYCACHE_H__
#define __SEQREGIONENTRYCACHE_H__

#include "EnsC.h"

typedef struct seqRegionCacheEntryStruct {
  char * regionName;
  IDType regionId;
  IDType csId;
  long   regionLength;
} SeqRegionCacheEntry;

SeqRegionCacheEntry *SeqRegionCacheEntry_new(IDType regionId, char *regionName, IDType csId, long regionLength);
void SeqRegionCacheEntry_free(SeqRegionCacheEntry *srce);

#endif
