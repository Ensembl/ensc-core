#include "Homology.h"
#include <stdio.h>

Homology *Homology_new() {
  Homology *hm;

  if ((hm = (Homology *)calloc(1,sizeof(Homology))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for hm\n");
    return NULL;
  }

  hm->objectType = CLASS_HOMOLOGY;
  return hm;
}

char *Homology_setSpecies(Homology *homol, char *species) {
  StrUtil_copyString(&(homol->species),species,0);
 
  return homol->species;
}

char *Homology_setStableId(Homology *homol, char *sid) {
  StrUtil_copyString(&(homol->stableId),sid,0);
 
  return homol->stableId;
}

char *Homology_setChromosome(Homology *homol, char *chr) {
  StrUtil_copyString(&(homol->chrName),chr,0);
 
  return homol->chrName;
}
