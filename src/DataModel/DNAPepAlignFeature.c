#include "DNAPepAlignFeature.h"

DNAPepAlignFeature *DNAPepAlignFeature_new() {
  DNAPepAlignFeature *dpaf;

  if ((dpaf = (DNAPepAlignFeature *)calloc(1,sizeof(DNAPepAlignFeature))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for dna pep align feature\n");
    return NULL;
  }

  dpaf->objectType = CLASS_DNAPEPALIGNFEATURE;

  dpaf->funcs = &baseAlignFeatureFuncs;
  return dpaf;
}
