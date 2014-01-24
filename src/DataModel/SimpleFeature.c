/*
 * Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

#define __SIMPLEFEATURE_MAIN__
#include "SimpleFeature.h"
#undef __SIMPLEFEATURE__MAIN__
#include "StrUtil.h"

#include <string.h>

SimpleFeature *SimpleFeature_new() {
  SimpleFeature *sf;

  if ((sf = (SimpleFeature *)calloc(1,sizeof(SimpleFeature))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for sf\n");
    return NULL;
  }

  sf->objectType = CLASS_SIMPLEFEATURE;
  Object_incRefCount(sf);

  sf->funcs = &simpleFeatureFuncs;
 
  return sf;
}

ECOSTRING SimpleFeature_setDisplayLabel(SimpleFeature *sf, char *label) {
  EcoString_copyStr(ecoSTable, &(sf->displayLabel), label, 0);
  return sf->displayLabel;
}

SimpleFeature *SimpleFeature_shallowCopyImpl(SimpleFeature *sf) {
  SimpleFeature *newSimpleFeature = SimpleFeature_new();

  memcpy(newSimpleFeature,sf,sizeof(SimpleFeature));

  return newSimpleFeature;
}

void SimpleFeature_freeImpl(SimpleFeature *sf) {
  Object_decRefCount(sf);

  if (Object_getRefCount(sf) > 0) {
    return;
  } else if (Object_getRefCount(sf) < 0) {
    fprintf(stderr,"Error: Negative reference count for SimpleFeature\n"
                   "       Freeing it anyway\n");
  }

  if (sf->displayLabel) EcoString_freeStr(ecoSTable, sf->displayLabel);

  SeqFeature_freePtrs((SeqFeature *)sf);
  free(sf);
}

