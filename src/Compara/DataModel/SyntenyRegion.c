/*
 * Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

#define __SYNTENYREGION_MAIN__
#include "SyntenyRegion.h"
#undef __SYNTENYREGION_MAIN__

SyntenyRegion *SyntenyRegion_new() {
  SyntenyRegion *sr;

  if ((sr = (SyntenyRegion *)calloc(1,sizeof(SyntenyRegion))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for sr\n");
    return NULL;
  }

  sr->objectType = CLASS_SYNTENYREGION;

  sr->funcs = &syntenyRegionFuncs;

  Object_incRefCount(sr);
  return sr;
}

ECOSTRING SyntenyRegion_setChrName(SyntenyRegion *sr, char *chrName) {
  EcoString_copyStr(ecoSTable, &(sr->chrName), chrName, 0);

  return sr->chrName;
}

ECOSTRING SyntenyRegion_setHitChrName(SyntenyRegion *sr, char *hitChrName) {
  EcoString_copyStr(ecoSTable, &(sr->hitChrName), hitChrName, 0);

  return sr->hitChrName;
}

ECOSTRING SyntenyRegion_setHitSeqType(SyntenyRegion *sr, char *hitSeqType) {
  EcoString_copyStr(ecoSTable, &(sr->hitSeqType), hitSeqType, 0);

  return sr->hitSeqType;
}

ECOSTRING SyntenyRegion_setSeqType(SyntenyRegion *sr, char *seqType) {
  EcoString_copyStr(ecoSTable, &(sr->seqType), seqType, 0);

  return sr->seqType;
}

void SyntenyRegion_free(SyntenyRegion *sr) {
  Object_decRefCount(sr);

  if (Object_getRefCount(sr) > 0) {
    return;
  } else if (Object_getRefCount(sr) < 0) {
    fprintf(stderr,"Error: Negative reference count for SyntenyRegion\n"
                   "       Freeing it anyway\n");
  }

  if (sr->chrName)    EcoString_freeStr(ecoSTable, sr->chrName);
  if (sr->hitChrName) EcoString_freeStr(ecoSTable, sr->hitChrName);
  if (sr->seqType)    EcoString_freeStr(ecoSTable, sr->seqType);
  if (sr->hitSeqType) EcoString_freeStr(ecoSTable, sr->hitSeqType);

  free(sr);
}

