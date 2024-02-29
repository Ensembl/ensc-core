/*
 * Copyright [1999-2024] EMBL-European Bioinformatics Institute
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

#define __REPEATCONSENSUS_MAIN__
#include "RepeatConsensus.h"
#undef __REPEATCONSENSUS_MAIN__
#include "EcoString.h"

RepeatConsensus *RepeatConsensus_new() {
  RepeatConsensus *rc;

  if ((rc = (RepeatConsensus *)calloc(1,sizeof(RepeatConsensus))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for rc\n");
    return NULL;
  }

  rc->objectType = CLASS_REPEATCONSENSUS;

  rc->funcs = &repeatConsensusFuncs;

  Object_incRefCount(rc);
  return rc;
}

ECOSTRING RepeatConsensus_setConsensus(RepeatConsensus *rc, char *cons) {
  EcoString_copyStr(ecoSTable,&(rc->consensus), cons, 0);
  return rc->consensus;
}

ECOSTRING RepeatConsensus_setRepeatClass(RepeatConsensus *rc, char *classStr) {
  EcoString_copyStr(ecoSTable,&(rc->repeatClass), classStr, 0);
  return rc->repeatClass;
}

ECOSTRING RepeatConsensus_setRepeatType(RepeatConsensus *rc, char *type) {
  EcoString_copyStr(ecoSTable,&(rc->repeatType), type, 0);
  return rc->repeatType;
}

ECOSTRING RepeatConsensus_setName(RepeatConsensus *rc, char *name) {
  EcoString_copyStr(ecoSTable,&(rc->name), name, 0);
  return rc->name;
}

void RepeatConsensus_free(RepeatConsensus *rc) {
  Object_decRefCount(rc);

  if (Object_getRefCount(rc) > 0) {
    return;
  } else if (Object_getRefCount(rc) < 0) {
    fprintf(stderr,"Error: Negative reference count for RepeatConsensus\n"
                   "       Freeing it anyway\n");
  }

  if (rc->name) EcoString_freeStr(ecoSTable, rc->name);
  if (rc->repeatClass) EcoString_freeStr(ecoSTable, rc->repeatClass);

  free(rc);
}
