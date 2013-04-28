#include "SeqRegionCacheEntry.h"
#include "StrUtil.h"

SeqRegionCacheEntry *SeqRegionCacheEntry_new(IDType regionId, char *regionName, IDType csId, long regionLength) {
  SeqRegionCacheEntry *cacheData;

  // Allocate the struct
  if ((cacheData = (SeqRegionCacheEntry *)calloc(1,sizeof(SeqRegionCacheEntry))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for SeqRegionCacheEntry\n");
    exit(1);
  }

  StrUtil_copyString(&(cacheData->regionName), regionName, 0);
  cacheData->regionId     = regionId;
  cacheData->csId         = csId;
  cacheData->regionLength = regionLength;

  return cacheData;
}

void SeqRegionCacheEntry_free(SeqRegionCacheEntry *srce) {
  free(srce->regionName);
  free(srce);
}
