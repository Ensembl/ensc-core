#define __FEATUREPAIR_MAIN__
#include "FeaturePair.h"
#undef __FEATUREPAIR_MAIN__
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

char *FeaturePair_setHitSeqName(FeaturePair *fp, char *str) {
  fp->hitId = StrUtil_copyString(&(fp->hitId),str,0);

  return fp->hitId;
}

char *FeaturePair_setHitSpecies(FeaturePair *fp, char *str) {
  fp->hitSpecies = StrUtil_copyString(&(fp->hitSpecies),str,0);

  return fp->hitSpecies;
}

char *FeaturePair_setSpecies(FeaturePair *fp, char *str) {
  fp->species = StrUtil_copyString(&(fp->species),str,0);

  return fp->species;
}

