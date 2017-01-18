/*
 * Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
 * Copyright [2016-2017] EMBL-European Bioinformatics Institute
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *      http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

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

