/*
 * Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

#ifndef __REPEATCONSENSUS_H__
#define __REPEATCONSENSUS_H__

#include "Storable.h"
#include "EcoString.h"
#include "EnsRoot.h"

OBJECTFUNC_TYPES(RepeatConsensus)

typedef struct RepeatConsensusFuncsStruct {
  OBJECTFUNCS_DATA(RepeatConsensus)
} RepeatConsensusFuncs;


#define FUNCSTRUCTTYPE RepeatConsensusFuncs
struct RepeatConsensusStruct {
  ENSROOT_DATA
  Storable st;
  ECOSTRING repeatClass;
  ECOSTRING consensus;
  ECOSTRING name;
  ECOSTRING repeatType;
  int length;
};
#undef FUNCSTRUCTTYPE

#define RepeatConsensus_setDbID(rcs,id) Storable_setDbID(&((rcs)->st),(id))
#define RepeatConsensus_getDbID(rcs) Storable_getDbID(&((rcs)->st))

#define RepeatConsensus_setAdaptor(rcs,ad) Storable_setAdaptor(&((rcs)->st),(ad))
#define RepeatConsensus_getAdaptor(rcs) Storable_getAdaptor(&((rcs)->st))

ECOSTRING RepeatConsensus_setRepeatClass(RepeatConsensus *rc, char *classStr);
#define RepeatConsensus_getRepeatClass(rcs) (rcs)->repeatClass

ECOSTRING RepeatConsensus_setRepeatType(RepeatConsensus *rc, char *type);
#define RepeatConsensus_getRepeatType(rcs) (rcs)->repeatType

ECOSTRING RepeatConsensus_setConsensus(RepeatConsensus *rc, char *cons);
#define RepeatConsensus_getConsensus(rcs) (rcs)->consensus

ECOSTRING RepeatConsensus_setName(RepeatConsensus *rc, char *name);
#define RepeatConsensus_getName(rcs) (rcs)->name

#define RepeatConsensus_setLength(rcs,len) (rcs)->length = (len)
#define RepeatConsensus_getLength(rcs) (rcs)->length

RepeatConsensus *RepeatConsensus_new();

void RepeatConsensus_free(RepeatConsensus *rc);

#ifdef __REPEATCONSENSUS_MAIN__
  RepeatConsensusFuncs
    repeatConsensusFuncs = {
                    RepeatConsensus_free,
                    NULL, // shallowCopy
                    NULL // deepCopy
                   };
#else
  extern RepeatConsensusFuncs repeatConsensusFuncs;
#endif


#endif
