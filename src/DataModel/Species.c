#include "Species.h"
#include "StrUtil.h"

Species *Species_new() {
  Species *species;

  if ((species = (Species *)calloc(1,sizeof(Species))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for species\n");
    return NULL;
  }

  return species;
}


char *Species_setCommonName(Species *species, char *comName) {
  StrUtil_copyString(&(species->commonName), comName,0);

  return species->commonName;
}

char *Species_getBinomialName(Species *species) {
  if (!species->binomialName) {
    char tmpStr[1024];
    if (Vector_getNumElement(species->classification) > 1) {
      sprintf(tmpStr, "%s %s\n", Vector_getElementAt(species->classification, 0),
                                 Vector_getElementAt(species->classification, 0));

    } else if (Vector_getNumElement(species->classification) == 1) {
      fprintf(stderr, "Warning: No genus info for binomial name\n");
      strcpy(tmpStr,Vector_getElementAt(species->classification,0));

    } else {
      strcpy(tmpStr,"sp.");
      fprintf(stderr, "Warning: No species info for binomial name\n");
    }
    StrUtil_copyString(&(species->binomialName), tmpStr, 0);
  }

  return species->binomialName;

}
