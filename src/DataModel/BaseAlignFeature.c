#include "BaseAlignFeature.h"

BaseAlignFeature *BaseAlignFeature_new() {
  BaseAlignFeature *baf;

  if ((baf = (BaseAlignFeature *)calloc(1,sizeof(BaseAlignFeature))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for base align feature\n");
    return NULL;
  }

  return baf;
}
