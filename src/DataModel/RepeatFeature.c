#include "RepeatFeature.h"

RepeatFeature *RepeatFeature_new() {
  RepeatFeature *rf;

  if ((rf = (RepeatFeature *)calloc(1,sizeof(RepeatFeature))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for repeat feature\n");
    return NULL;
  }

  rf->objectType = CLASS_REPEATFEATURE;

  rf->funcs = &seqFeatureFuncs;

  return rf;
}
