#include "GenomicRange.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

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

void GenomicRange_expand(GenomicRange *range, int pad) {
  GenomicRange_setChrStart(range,GenomicRange_getChrStart(range) - pad);
  GenomicRange_setChrEnd(range,GenomicRange_getChrEnd(range) + pad);
  
  if (GenomicRange_getChrStart(range) < 1) {
    GenomicRange_setChrStart(range,1); 
  } 
}

void GenomicRange_free(GenomicRange *range) {
  free(range->chrName);
  free(range);
}
