#include "Chromosome.h"

Chromosome *Chromosome_new() {
  Chromosome *chrom;

  if ((chrom = (Chromosome *)calloc(1,sizeof(Chromosome))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for chrom\n");
    return NULL;
  }

  return chrom;
}
