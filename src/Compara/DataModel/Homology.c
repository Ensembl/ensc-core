/*
 * Copyright [1999-2024] EMBL-European Bioinformatics Institute
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

#define __HOMOLOGY_MAIN__
#include "Homology.h"
#undef __HOMOLOGY_MAIN__
#include <stdio.h>

Homology *Homology_new() {
  Homology *hm;

  if ((hm = (Homology *)calloc(1,sizeof(Homology))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for hm\n");
    return NULL;
  }
  Object_incRefCount(hm);

  hm->objectType = CLASS_HOMOLOGY;

  hm->funcs = &homologyFuncs;

  return hm;
}

ECOSTRING Homology_setSpecies(Homology *homol, char *species) {
  EcoString_copyStr(ecoSTable, &(homol->species),species,0);
 
  return homol->species;
}

ECOSTRING Homology_setStableId(Homology *homol, char *sid) {
  EcoString_copyStr(ecoSTable, &(homol->stableId),sid,0);
 
  return homol->stableId;
}

ECOSTRING Homology_setChromosome(Homology *homol, char *chr) {
  EcoString_copyStr(ecoSTable, &(homol->chrName),chr,0);
 
  return homol->chrName;
}

void Homology_free(Homology *hom) {
  Object_decRefCount(hom);

  if (Object_getRefCount(hom) > 0) {
    return;
  } else if (Object_getRefCount(hom) < 0) {
    fprintf(stderr,"Error: Negative reference count for Homology\n"
                   "       Freeing it anyway\n");
  }

  if (hom->chrName)  EcoString_freeStr(ecoSTable, hom->chrName);
  if (hom->stableId) EcoString_freeStr(ecoSTable, hom->stableId);
  if (hom->species)  EcoString_freeStr(ecoSTable, hom->species);

  free(hom);
}

