#include "GenomicAlign.h"

GenomicAlign *GenomicAlign_new() {
  GenomicAlign *ga;

  if ((ga = (GenomicAlign *)calloc(1,sizeof(GenomicAlign))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for ga\n");
    return NULL;
  }

  ga->objectType = CLASS_GENOMICALIGN;
  return ga;
}

