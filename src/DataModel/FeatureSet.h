/*
 * Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
 * Copyright [2016-2017] EMBL-European Bioinformatics Institute
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *      http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

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
void FeatureSet_removeAll(FeatureSet *fs);


#endif
