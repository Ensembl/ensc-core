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

#define __ATTRIBUTE_MAIN__
#include "Attribute.h"
#undef __ATTRIBUTE_MAIN__
#include "StrUtil.h"
#include <string.h>
#include "EcoString.h"

Attribute *Attribute_new() {
  Attribute *attribute;

  if ((attribute = (Attribute *)calloc(1,sizeof(Attribute))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for attribute\n");
    return NULL;
  }

  attribute->funcs = &attributeFuncs;

  attribute->objectType = CLASS_ATTRIBUTE;

  return attribute;
}

char *Attribute_setValue(Attribute *attribute, char *value) {
  StrUtil_copyString(&(attribute->value), value,0);

  return attribute->value;
}

ECOSTRING Attribute_setCode(Attribute *attrib, char *code) {
  EcoString_copyStr(ecoSTable, &(attrib->code),code,0);

  if (attrib->code == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for code\n");
    return NULL;
  }

  return attrib->code;
}

ECOSTRING Attribute_setName(Attribute *attrib, char *name) {
  EcoString_copyStr(ecoSTable, &(attrib->name),name,0);

  if (attrib->name == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for name\n");
    return NULL;
  }

  return attrib->name;
}

ECOSTRING Attribute_setDescription(Attribute *attrib, char *description) {
  EcoString_copyStr(ecoSTable, &(attrib->description),description,0);

  if (attrib->description == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for description\n");
    return NULL;
  }

  return attrib->description;
}

void Attribute_free(Attribute *attribute) {
  EcoString_freeStr(ecoSTable, attribute->name);
  EcoString_freeStr(ecoSTable, attribute->description);
  EcoString_freeStr(ecoSTable, attribute->code);

  free(attribute->value);

  free(attribute);
}
