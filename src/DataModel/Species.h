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
                    Species_free
                   };
#else
  extern SpeciesFuncs speciesFuncs;
#endif



#endif
