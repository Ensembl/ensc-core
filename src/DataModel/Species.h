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

#ifndef __SPECIES_H__
#define __SPECIES_H__

#include "DataModelTypes.h"
#include "Vector.h"
#include "EnsRoot.h"

OBJECTFUNC_TYPES(Species)

typedef struct SpeciesFuncsStruct {
  OBJECTFUNCS_DATA(Species)
} SpeciesFuncs;

#define FUNCSTRUCTTYPE SpeciesFuncs
struct SpeciesStruct {
  ENSROOT_DATA
  Vector *classification;
  char *commonName;
  char *binomialName;
};
#undef FUNCSTRUCTTYPE

Species *Species_new();

#define Species_setClassification(sp, cl) (sp)->classification = (cl)
#define Species_getClassification(sp) (sp)->classification

char *Species_setCommonName(Species *species, char *commonName);
#define Species_getCommonName(sp) (sp)->commonName

char *Species_getBinomialName(Species *species);

void Species_free(Species *species);

#ifdef __SPECIES_MAIN__
  SpeciesFuncs
    speciesFuncs = {
                    Species_free,
                    NULL, // shallowCopy
                    NULL  // deepCopy
                   };
#else
  extern SpeciesFuncs speciesFuncs;
#endif



#endif
