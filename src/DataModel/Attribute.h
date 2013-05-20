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
