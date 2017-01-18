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

#define __DNAALIGNFEATURE_MAIN__
#include "DNAAlignFeature.h"
#undef __DNAALIGNFEATURE_MAIN__
//#include "ProcUtil.h"

#include <string.h>

DNAAlignFeature *DNAAlignFeature_new() {
  DNAAlignFeature *daf;

  if ((daf = (DNAAlignFeature *)calloc(1,sizeof(DNAAlignFeature))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for dna align feature\n");
    return NULL;
  }

  daf->objectType = CLASS_DNADNAALIGNFEATURE;
//  Object_incRefCount(daf);

  daf->funcs = &dnaAlignFeatureFuncs;

// Not very happy with this way of signifying values not set, but for now can't think of a better way which is efficient
  DNAAlignFeature_setpValue(daf, FLOAT_UNDEF);
  DNAAlignFeature_setPercId(daf, FLOAT_UNDEF);
  DNAAlignFeature_sethCoverage(daf, FLOAT_UNDEF);

  return daf;
}

int DNAAlignFeature_getHitUnit(void) {
  return 1;
}

int DNAAlignFeature_getQueryUnit(void) {
  return 1;
}

ECOSTRING DNAAlignFeature_setExtraData(DNAAlignFeature *daf, char *extraData) {
  EcoString_copyStr(ecoSTable, (char**)(&(daf->extraData)), extraData,0);

  if (daf->extraData == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for extraData\n");
    return NULL;
  }

  return daf->extraData;
}

void DNAAlignFeature_freeImpl(DNAAlignFeature *daf) {
  Object_decRefCount(daf);
//  Object_decRefCount(daf);

//  fprintf(stderr,"DAF_freeImpl called ref count %d\n", Object_decRefCount(daf));
  if (Object_getRefCount(daf) > 0) {
//    ProcUtil_showBacktrace(EnsC_progName);
//    printf("return\n");
    return;
  } else if (Object_getRefCount(daf) < 0) {
//    fprintf(stderr,"Error: Negative reference count for DNAAlignFeature\n"
//                   "       Freeing it anyway\n");
  }

  BaseAlignFeature_freePtrs((BaseAlignFeature *)daf);
//  fprintf(stderr," freeing daf\n");
  
  free(daf);
}

DNAAlignFeature *DNAAlignFeature_shallowCopyImpl(DNAAlignFeature *daf) {
  DNAAlignFeature *newDNAAlignFeature = DNAAlignFeature_new();

  memcpy(newDNAAlignFeature,daf,sizeof(DNAAlignFeature));

  return newDNAAlignFeature;
}

DNAAlignFeature *DNAAlignFeature_deepCopyImpl(DNAAlignFeature *daf) {
  DNAAlignFeature *newDNAAlignFeature = DNAAlignFeature_new();

  memcpy(newDNAAlignFeature,daf,sizeof(DNAAlignFeature));

  BaseAlignFeature_copyData((BaseAlignFeature*)newDNAAlignFeature, (BaseAlignFeature*)daf);
  

  newDNAAlignFeature->referenceCount = 0;

  return newDNAAlignFeature;
}
