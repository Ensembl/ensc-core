#define __BASEALIGNFEATURE_MAIN__
#include "FeaturePair.h"
#undef __BASEALIGNFEATURE_MAIN__
#include "StrUtil.h"

FeaturePair *FeaturePair_new(void) {
  FeaturePair *fp;

  if ((fp = (FeaturePair *)calloc(1,sizeof(FeaturePair))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for feature pair\n");
    return NULL;
  }

  fp->objectType = CLASS_FEATUREPAIR;

  fp->funcs = &featurePairFuncs;

  return fp;
}

char *FeaturePair_setHitId(FeaturePair *fp, char *str) {
  fp->hitId = StrUtil_copyString(&(fp->hitId),str,0);

  return fp->hitId;
}

