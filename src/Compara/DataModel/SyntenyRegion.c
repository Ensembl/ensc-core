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

