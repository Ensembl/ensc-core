#define __ATTRIBUTE_MAIN__
#include "Attribute.h"
#undef __ATTRIBUTE_MAIN__
#include "StrUtil.h"
#include <string.h>

Attribute *Attribute_new() {
  Attribute *attribute;

  if ((attribute = (Attribute *)calloc(1,sizeof(Attribute))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for attribute\n");
    return NULL;
  }

  attribute->funcs = &attributeFuncs;

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
  fprintf(stderr,"Warning: Attribute_free NIY\n");
}
