#ifndef __FEATURESET_H__
#define __FEATURESET_H__

#include "DataModelTypes.h"
#include "EnsC.h"

struct FeatureSetStruct {
  void **features;
  int nFeature;
};

void *FeatureSet_addFeature(FeatureSet *fset, void *sf);
#define FeatureSet_getNumFeature(fs) (fs)->nFeature
#define FeatureSet_getFeatures(fs) (fs)->features
void *FeatureSet_getFeatureAt(FeatureSet *fset,int ind);

#endif
