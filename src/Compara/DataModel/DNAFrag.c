#include "DNAFrag.h"

DNAFrag *DNAFrag_new() {
  DNAFrag *df;

  if ((df = (DNAFrag *)calloc(1,sizeof(DNAFrag))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for df\n");
    return NULL;
  }

  df->objectType = CLASS_DNAFRAG;
  return df;
}

