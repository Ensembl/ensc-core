#include "Homology.h"

Homology *Homology_new() {
  Homology *hm;

  if ((hm = (Homology *)calloc(1,sizeof(Homology))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for hm\n");
    return NULL;
  }

  hm->objectType = CLASS_HOMOLOGY;
  return hm;
}

