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

#ifndef __STABLEIDINFO_H__
#define __STABLEIDINFO_H__

#include <time.h>
#include <stdio.h>

#include "DataModelTypes.h"

struct StableIdInfoStruct {
  char *stableId;
  int   version;
  time_t created;
  time_t modified;
  char isCurrent;
};

char *StableIdInfo_setStableId(StableIdInfo *si, char *sid); 
#define StableIdInfo_getStableId(si)  (si)->stableId

#define StableIdInfo_setCreated(si,cd)  (si)->created = (cd)
#define StableIdInfo_getCreated(si)  (si)->created

#define StableIdInfo_setModified(si,mod)  (si)->modified = (mod)
#define StableIdInfo_getModified(si) (si)->modified

#define StableIdInfo_setVersion(si,ver)  (si)->version = (ver)
#define StableIdInfo_getVersion(si) (si)->version

#define StableIdInfo_setIsCurrent(si,isC)  (si)->isCurrent = (isC)
#define StableIdInfo_getIsCurrent(si) (si)->isCurrent

void StableIdInfo_freePtrs(StableIdInfo *si);

#endif
