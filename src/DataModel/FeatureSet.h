#ifndef __FEATURESET_H__
#define __FEATURESET_H__

#include "EnsC.h"

typedef struct FeatureSetStruct {
  void **features;
  int nFeature;
} FeatureSet;

void *FeatureSet_addFeature(FeatureSet *fset, void *sf);
#define FeatureSet_getNumFeature(fs) (fs)->nFeature
void *FeatureSet_getFeatureAt(FeatureSet *fset,int ind);

#endif
