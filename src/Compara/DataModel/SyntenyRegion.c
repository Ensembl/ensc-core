#define __SYNTENYREGION_MAIN__
#include "SyntenyRegion.h"
#undef __SYNTENYREGION_MAIN__

SyntenyRegion *SyntenyRegion_new() {
  SyntenyRegion *sr;

  if ((sr = (SyntenyRegion *)calloc(1,sizeof(SyntenyRegion))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for sr\n");
    return NULL;
  }

  sr->objectType = CLASS_SYNTENYREGION;

  sr->funcs = &syntenyRegionFuncs;

  Object_incRefCount(sr);
  return sr;
}

ECOSTRING SyntenyRegion_setChrName(SyntenyRegion *sr, char *chrName) {
  EcoString_copyStr(ecoSTable, &(sr->chrName), chrName, 0);

  return sr->chrName;
}

ECOSTRING SyntenyRegion_setHitChrName(SyntenyRegion *sr, char *hitChrName) {
  EcoString_copyStr(ecoSTable, &(sr->hitChrName), hitChrName, 0);

  return sr->hitChrName;
}

ECOSTRING SyntenyRegion_setHitSeqType(SyntenyRegion *sr, char *hitSeqType) {
  EcoString_copyStr(ecoSTable, &(sr->hitSeqType), hitSeqType, 0);

  return sr->hitSeqType;
}

ECOSTRING SyntenyRegion_setSeqType(SyntenyRegion *sr, char *seqType) {
  EcoString_copyStr(ecoSTable, &(sr->seqType), seqType, 0);

  return sr->seqType;
}

void SyntenyRegion_free(SyntenyRegion *sr) {
  Object_decRefCount(sr);

  if (Object_getRefCount(sr) > 0) {
    return;
  } else if (Object_getRefCount(sr) < 0) {
    fprintf(stderr,"Error: Negative reference count for SyntenyRegion\n"
                   "       Freeing it anyway\n");
  }

  if (sr->chrName)    EcoString_freeStr(ecoSTable, sr->chrName);
  if (sr->hitChrName) EcoString_freeStr(ecoSTable, sr->hitChrName);
  if (sr->seqType)    EcoString_freeStr(ecoSTable, sr->seqType);
  if (sr->hitSeqType) EcoString_freeStr(ecoSTable, sr->hitSeqType);

  free(sr);
}

