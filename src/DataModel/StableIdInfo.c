#include "StableIdInfo.h"

#include <stdlib.h>
#include <string.h>

char *StableIdInfo_setStableId(StableIdInfo *si, char *sid) {
  if ((si->stableId = (char *)malloc(strlen(sid)+1)) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for stableid\n");
    return NULL;
  }

  strcpy(si->stableId,sid);

  return si->stableId;
}

void StableIdInfo_freePtrs(StableIdInfo *si) {
  if (si->stableId) free(si->stableId);
}
