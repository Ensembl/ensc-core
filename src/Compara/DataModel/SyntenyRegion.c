#include "SyntenyRegion.h"

SyntenyRegion *SyntenyRegion_new() {
  SyntenyRegion *sr;

  if ((sr = (SyntenyRegion *)calloc(1,sizeof(SyntenyRegion))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for sr\n");
    return NULL;
  }

  sr->objectType = CLASS_SYNTENYREGION;
  return sr;
}

char *SyntenyRegion_setChrName(SyntenyRegion *sr, char *chrName) {
  StrUtil_copyString(&(sr->chrName), chrName, 0);

  return sr->chrName;
}

char *SyntenyRegion_setHitChrName(SyntenyRegion *sr, char *hitChrName) {
  StrUtil_copyString(&(sr->hitChrName), hitChrName, 0);

  return sr->hitChrName;
}

char *SyntenyRegion_setHitSeqType(SyntenyRegion *sr, char *hitSeqType) {
  StrUtil_copyString(&(sr->hitSeqType), hitSeqType, 0);

  return sr->hitSeqType;
}

char *SyntenyRegion_setSeqType(SyntenyRegion *sr, char *seqType) {
  StrUtil_copyString(&(sr->seqType), seqType, 0);

  return sr->seqType;
}
