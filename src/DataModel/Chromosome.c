#define __CHROMOSOME_MAIN__
#include "Chromosome.h"
#undef __CHROMOSOME_MAIN__

Chromosome *Chromosome_new() {
  Chromosome *chrom;

  if ((chrom = (Chromosome *)calloc(1,sizeof(Chromosome))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for chrom\n");
    return NULL;
  }

  chrom->objectType = CLASS_CHROMOSOME;

  chrom->funcs = &chromosomeFuncs;

  return chrom;
}

ECOSTRING Chromosome_setName(Chromosome *c,char *name) {
  EcoString_copyStr(ecoSTable, &(c->name), name, 0);

  return c->name;
}

void Chromosome_free(Chromosome *chr) {
  Object_decRefCount(chr);

  if (Object_getRefCount(chr) > 0) {
    return;
  } else if (Object_getRefCount(chr) < 0) {
    fprintf(stderr,"Error: Negative reference count for Chromosome\n"
                   "       Freeing it anyway\n");
  }

  if (chr->name) EcoString_freeStr(ecoSTable, chr->name);

  free(chr);
}

