#define __REPEATFEATURE_MAIN__
#include "RepeatFeature.h"
#undef __REAPEATFEATURE_MAIN__

RepeatFeature *RepeatFeature_new() {
  RepeatFeature *rf;

  if ((rf = (RepeatFeature *)calloc(1,sizeof(RepeatFeature))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for repeat feature\n");
    return NULL;
  }

  rf->objectType = CLASS_REPEATFEATURE;
  Object_incRefCount(rf);

  rf->funcs = &repeatFeatureFuncs;

  return rf;
}

void RepeatFeature_free(RepeatFeature *rf) {
  Object_decRefCount(rf);

  if (Object_getRefCount(rf) > 0) {
    return;
  } else if (Object_getRefCount(rf) < 0) {
    fprintf(stderr,"Error: Negative reference count for RepeatFeature\n"
                   "       Freeing it anyway\n");
  }

  SeqFeature_freePtrs((SeqFeature *)rf);
  free(rf);
}

