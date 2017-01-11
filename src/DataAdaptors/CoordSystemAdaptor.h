/*
 * Copyright [1999-2017] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

#ifndef __COORDSYSTEMADAPTOR_H__
#define __COORDSYSTEMADAPTOR_H__

#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "CoordSystem.h"
#include "IDHash.h"
#include "StringHash.h"
#include "Vector.h"

struct CoordSystemAdaptorStruct {
  BASEADAPTOR_DATA
  IDHash *dbIDCache;
  IDHash *rankCache;
  StringHash *nameCache;

  IDHash *isSeqLevelCache;
  IDHash *isDefaultVersionCache;

  StringHash *mappingPaths;

  CoordSystem *topLevel;
};


int CoordSystem_sortByRankFunc(const void *one, const void *two);

void                CoordSystemAdaptor_cacheMappingPaths(CoordSystemAdaptor *csa);
void                CoordSystemAdaptor_cacheSeqRegionMapping(CoordSystemAdaptor *csa);
void                CoordSystemAdaptor_dumpCachedMappings(CoordSystemAdaptor *csa);
Vector *            CoordSystemAdaptor_fetchAllByAttrib(CoordSystemAdaptor *csa, char *attrib);
Vector *            CoordSystemAdaptor_fetchAllByName(CoordSystemAdaptor *csa, char *name);
CoordSystem *       CoordSystemAdaptor_fetchByAttrib(CoordSystemAdaptor *csa, char *attrib, char *version);
CoordSystem *       CoordSystemAdaptor_fetchByDbID(CoordSystemAdaptor *csa, IDType dbID);
CoordSystem *       CoordSystemAdaptor_fetchByName(CoordSystemAdaptor *csa, char *name, char *version);
CoordSystem *       CoordSystemAdaptor_fetchByRank(CoordSystemAdaptor *csa, int rank);
Vector *            CoordSystemAdaptor_fetchAll(CoordSystemAdaptor *csa);
CoordSystem *       CoordSystemAdaptor_fetchSeqLevel(CoordSystemAdaptor *csa);
CoordSystem *       CoordSystemAdaptor_fetchTopLevel(CoordSystemAdaptor *csa);
Vector *            CoordSystemAdaptor_getMappingPath(CoordSystemAdaptor *csa, CoordSystem *cs1, CoordSystem *cs2);
CoordSystemAdaptor *CoordSystemAdaptor_new(DBAdaptor *dba);

#endif
