#include "Species.h"

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

  return NULL;
}
