#include "FeatureSet.h"

void *FeatureSet_getFeatureAt(FeatureSet *fs, int ind) {
  if (ind < 0 || ind >= fs->nFeature) {
    fprintf(stderr,"ERROR: Invalid feature index %d\n",ind);
    exit(1);
  }
  return fs->features[ind];
}

void *FeatureSet_addFeature(FeatureSet *fs, void *sf) {
  if (sf == NULL) {
    fprintf(stderr, "ERROR: Feature null in FeatureSet_addFeature call\n");
    return NULL;
  }

  if (!fs->nFeature) fs->features = NULL;

  fs->nFeature++;
  if ((fs->features = (void **)realloc(fs->features,
               fs->nFeature*sizeof(void *))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for sf array\n");
    return NULL;
  }

  fs->features[fs->nFeature-1] = sf;

  return sf;
}

