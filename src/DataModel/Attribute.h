/*
 * Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

#ifndef __ATTRIBUTE_H__
#define __ATTRIBUTE_H__

#include "DataModelTypes.h"
#include "Vector.h"
#include "EnsRoot.h"
#include "EcoString.h"

OBJECTFUNC_TYPES(Attribute)

typedef struct AttributeFuncsStruct {
  OBJECTFUNCS_DATA(Attribute)
} AttributeFuncs;

#define FUNCSTRUCTTYPE AttributeFuncs
struct AttributeStruct {
  ENSROOT_DATA
  char *value;
  ECOSTRING name;
  ECOSTRING code;
  ECOSTRING description;
};
#undef FUNCSTRUCTTYPE

Attribute *Attribute_new();

ECOSTRING Attribute_setName(Attribute *attribute, char *name);
#define Attribute_getName(attrib) (attrib)->name

ECOSTRING Attribute_setCode(Attribute *attribute, char *code);
#define Attribute_getCode(attrib) (attrib)->code

ECOSTRING Attribute_setDescription(Attribute *attribute, char *description);
#define Attribute_getDescription(attrib) (attrib)->description

ECOSTRING Attribute_setValue(Attribute *attribute, char *value);
#define Attribute_getValue(attrib) (attrib)->value

void Attribute_free(Attribute *attribute);

#ifdef __ATTRIBUTE_MAIN__
  AttributeFuncs
    attributeFuncs = {
                    Attribute_free,
                    NULL, // shallowCopy
                    NULL  // deepCopy
                   };
#else
  extern AttributeFuncs attributeFuncs;
#endif



#endif
