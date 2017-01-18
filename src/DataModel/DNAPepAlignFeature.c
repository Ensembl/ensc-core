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

#define __DNAPEPALIGNFEATURE_MAIN__
#include "DNAPepAlignFeature.h"
#undef __DNAPEPALIGNFEATURE_MAIN__

#include <string.h>

DNAPepAlignFeature *DNAPepAlignFeature_new() {
  DNAPepAlignFeature *dpaf;

  if ((dpaf = (DNAPepAlignFeature *)calloc(1,sizeof(DNAPepAlignFeature))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for dna pep align feature\n");
    return NULL;
  }

  dpaf->objectType = CLASS_DNAPEPALIGNFEATURE;
  Object_incRefCount(dpaf);

  dpaf->funcs = &dnaPepAlignFeatureFuncs;
  return dpaf;
}

int DNAPepAlignFeature_getHitUnit(void) {
  return 1; 
}

int DNAPepAlignFeature_getQueryUnit(void) {
  return 3;
}

void DNAPepAlignFeature_freeImpl(DNAPepAlignFeature *dpaf) {
  Object_decRefCount(dpaf);

  if (Object_getRefCount(dpaf) > 0) {
    return;
  } else if (Object_getRefCount(dpaf) < 0) {
    fprintf(stderr,"Error: Negative reference count for DNAPepAlignFeature\n"
                   "       Freeing it anyway\n");
  }

  BaseAlignFeature_freePtrs((BaseAlignFeature *)dpaf);

  free(dpaf);
}

DNAPepAlignFeature *DNAPepAlignFeature_shallowCopyImpl(DNAPepAlignFeature *dpaf) {
  DNAPepAlignFeature *newDNAPepAlignFeature = DNAPepAlignFeature_new();

  memcpy(newDNAPepAlignFeature,dpaf,sizeof(DNAPepAlignFeature));

  return newDNAPepAlignFeature;
}

DNAPepAlignFeature *DNAPepAlignFeature_deepCopyImpl(DNAPepAlignFeature *dpaf) {
  DNAPepAlignFeature *newDNAPepAlignFeature = DNAPepAlignFeature_new();

  memcpy(newDNAPepAlignFeature,dpaf,sizeof(DNAPepAlignFeature));

  return newDNAPepAlignFeature;
}
