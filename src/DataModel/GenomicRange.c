#include "GenomicRange.h"
#include <stdio.h>
#include <string.h>

GenomicRange *GenomicRange_new() {
  GenomicRange *genRange;

  if ((genRange = (GenomicRange *)calloc(1,sizeof(GenomicRange))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for GenomicRange\n");
    return NULL;
  }

  return genRange;
}

char *GenomicRange_setChrName(GenomicRange *gr, char *chrName) {
  if ((gr->chrName = (char *)malloc(strlen(chrName)+1)) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for chrName\n");
    exit(1);
  }

  strcpy(gr->chrName,chrName);

  return gr->chrName;
}

void GenomicRange_free(GenomicRange *range) {
  free(range->chrName);
  free(range);
}
