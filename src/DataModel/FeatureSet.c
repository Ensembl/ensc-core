/*
 * Copyright [1999-2017] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

void FeatureSet_removeAll(FeatureSet *fs) {
  fs->nFeature = 0;
}
