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

#define __FEATUREPAIR_MAIN__
#include "FeaturePair.h"
#undef __FEATUREPAIR_MAIN__
#include "StrUtil.h"

FeaturePair *FeaturePair_new(void) {
  FeaturePair *fp;

  if ((fp = (FeaturePair *)calloc(1,sizeof(FeaturePair))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for feature pair\n");
    return NULL;
  }

  fp->objectType = CLASS_FEATUREPAIR;
  Object_incRefCount(fp);

  fp->funcs = &featurePairFuncs;

  return fp;
}

void FeaturePair_copyData(FeaturePair *to, FeaturePair *from) {
  if (FeaturePair_getHitSeqName(from))
    FeaturePair_setHitSeqName(to, FeaturePair_getHitSeqName(from));
}

ECOSTRING FeaturePair_setHitSeqName(FeaturePair *fp, char *str) {
  EcoString_copyStr(ecoSTable,&(fp->hitId),str,0);

  return fp->hitId;
}

ECOSTRING FeaturePair_setHitSpecies(FeaturePair *fp, char *str) {
  EcoString_copyStr(ecoSTable, &(fp->hitSpecies),str,0);

  return fp->hitSpecies;
}

ECOSTRING FeaturePair_setSpecies(FeaturePair *fp, char *str) {
  EcoString_copyStr(ecoSTable,&(fp->species),str,0);

  return fp->species;
}

void FeaturePair_freePtrs(FeaturePair *fp) {
  if (fp->species) EcoString_freeStr(ecoSTable, fp->species);
  if (fp->hitSpecies) EcoString_freeStr(ecoSTable, fp->hitSpecies);
  if (fp->hitId) EcoString_freeStr(ecoSTable, fp->hitId);

  SeqFeature_freePtrs((SeqFeature *)fp);
}


void FeaturePair_freeImpl(FeaturePair *fp) {
  Object_decRefCount(fp);

  if (Object_getRefCount(fp) > 0) {
    return;
  } else if (Object_getRefCount(fp) < 0) {
    fprintf(stderr,"Error: Negative reference count for FeaturePair\n"
                   "       Freeing it anyway\n");
  }

  FeaturePair_freePtrs(fp);

  free(fp);
}
