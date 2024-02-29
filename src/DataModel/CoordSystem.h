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

#ifndef __COORDSYSTEM_H__
#define __COORDSYSTEM_H__

#include "DataModelTypes.h"
#include "Storable.h"

#include "EnsRoot.h"

OBJECTFUNC_TYPES(CoordSystem)

typedef struct CoordSystemFuncsStruct {
  OBJECTFUNCS_DATA(CoordSystem)
} CoordSystemFuncs;

#define FUNCSTRUCTTYPE CoordSystemFuncs
struct CoordSystemStruct {
  ENSROOT_DATA
  Storable st;
  char *name;
  char *version;
  char *nameColonVersion;
  char *dbIdStr;
  int   lenNameColonVersion;
  int   rank;

  int  isTopLevel;
  int  isSeqLevel;
  int  isDefaultVersion;
};
#undef FUNCSTRUCTTYPE

CoordSystem *CoordSystem_new(char *name, char *version, int rank, IDType dbID, CoordSystemAdaptor *csa, int isDefault, int isSeqLevel, int isTopLevel);

#define CoordSystem_setDbID(cs,dbID) Storable_setDbID(&((cs)->st),dbID)
#define CoordSystem_getDbID(cs) Storable_getDbID(&((cs)->st))

#define CoordSystem_setAdaptor(cs,ad) Storable_setAdaptor(&((cs)->st),ad)
#define CoordSystem_getAdaptor(cs) Storable_getAdaptor(&((cs)->st))

#define CoordSystem_getName(cs) cs->name

#define CoordSystem_getRank(cs) cs->rank

#define CoordSystem_getVersion(cs) cs->version

#define CoordSystem_getIsTopLevel(cs) cs->isTopLevel

#define CoordSystem_getIsSeqLevel(cs) cs->isSeqLevel

#define CoordSystem_getIsDefaultVersion(cs) cs->isDefaultVersion

int CoordSystem_compare(CoordSystem *cs1, CoordSystem *cs2);

char *CoordSystem_getNameColonVersion(CoordSystem *cs);
int CoordSystem_getLenNameColonVersion(CoordSystem *cs);

char *CoordSystem_getDbIDStr(CoordSystem *cs);

void CoordSystem_free(CoordSystem *coordSystem);

#ifdef __COORDSYSTEM_MAIN__
  CoordSystemFuncs
    coordSystemFuncs = {
                        CoordSystem_free,
                        NULL, // shallowCopy
                        NULL  // deepCopy
                       };
#else
  extern CoordSystemFuncs coordSystemFuncs;
#endif


#endif
