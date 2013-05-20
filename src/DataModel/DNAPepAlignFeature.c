#define __DNAPEPALIGNFEATURE_MAIN__
#include "DNAPepAlignFeature.h"
#undef __DNAPEPALIGNFEATURE_MAIN__

DNAPepAlignFeature *DNAPepAlignFeature_new() {
  DNAPepAlignFeature *dpaf;

  if ((dpaf = (DNAPepAlignFeature *)calloc(1,sizeof(DNAPepAlignFeature))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for dna pep align feature\n");
    return NULL;
  }

  dpaf->objectType = CLASS_DNAPEPALIGNFEATURE;
  Object_incRefCount(dpaf);

  dpaf->funcs = &dnaPepAlignFeatureFuncs;
  return dpaf;
}

int DNAPepAlignFeature_getHitUnit(void) {
  return 1; 
}

int DNAPepAlignFeature_getQueryUnit(void) {
  return 3;
}

void DNAPepAlignFeature_freeImpl(DNAPepAlignFeature *dpaf) {
  Object_decRefCount(dpaf);

  if (Object_getRefCount(dpaf) > 0) {
    return;
  } else if (Object_getRefCount(dpaf) < 0) {
    fprintf(stderr,"Error: Negative reference count for DNAPepAlignFeature\n"
                   "       Freeing it anyway\n");
  }

  BaseAlignFeature_freePtrs((BaseAlignFeature *)dpaf);

  free(dpaf);
}

DNAPepAlignFeature *DNAPepAlignFeature_shallowCopyImpl(DNAPepAlignFeature *dpaf) {
  DNAPepAlignFeature *newDNAPepAlignFeature = DNAPepAlignFeature_new();

  memcpy(newDNAPepAlignFeature,dpaf,sizeof(DNAPepAlignFeature));

  return newDNAPepAlignFeature;
}
