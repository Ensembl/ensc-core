#include "BaseAlignFeature.h"
#include "StrUtil.h"

BaseAlignFeature *BaseAlignFeature_new() {
  BaseAlignFeature *baf;

  if ((baf = (BaseAlignFeature *)calloc(1,sizeof(BaseAlignFeature))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for base align feature\n");
    return NULL;
  }

  return baf;
}

char *BaseAlignFeature_setCigarString(BaseAlignFeature *baf, char *str) {
  baf->cigarString = StrUtil_copyString(&(baf->cigarString),str,0);

  return baf->cigarString;
}

char *BaseAlignFeature_setHitId(BaseAlignFeature *baf, char *str) {
  baf->hitId = StrUtil_copyString(&(baf->cigarString),str,0);

  return baf->hitId;
}
