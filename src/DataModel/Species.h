#ifndef __SPECIES_H__
#define __SPECIES_H__

#include "DataModelTypes.h"
#include "Vector.h"
#include "EnsRoot.h"

#define FUNCSTRUCTTYPE NoTypeFuncs
struct SpeciesStruct {
  ENSROOT_DATA
  Vector *classification;
  char *commonName;
};
#undef FUNCSTRUCTTYPE

Species *Species_new();

#define Species_setClassification(sp, cl) (sp)->classification = (cl)
#define Species_getClassification(sp) (sp)->classification

char *Species_setCommonName(Species *species, char *commonName);
#define Species_getCommonName(sp) (sp)->commonName

char *Species_getBinomialName(Species *species);


#endif
