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

char *GenomicAlign_setCigarString(GenomicAlign *ga, char *cigStr) {
  StrUtil_copyString(&(ga->cigarString), cigStr, 0);

  return ga->cigarString;
}

char *GenomicAlign_setAlignmentType(GenomicAlign *ga, char *alType) {
  StrUtil_copyString(&(ga->alignmentType), alType, 0);

  return ga->alignmentType;
}
