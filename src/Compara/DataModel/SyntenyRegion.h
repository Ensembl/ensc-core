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

#ifndef __SYNTENYREGION_H__
#define __SYNTENYREGION_H__

#include "ComparaDataModelTypes.h"
#include "EnsRoot.h"
#include "Storable.h"

OBJECTFUNC_TYPES(SyntenyRegion)

typedef struct SyntenyRegionFuncsStruct {
  OBJECTFUNCS_DATA(SyntenyRegion)
} SyntenyRegionFuncs;

#define FUNCSTRUCTTYPE SyntenyRegionFuncs
struct SyntenyRegionStruct {
  ENSROOT_DATA
  Storable st;
  int start;
  int end;
  IDType clusterId;
  IDType dnaFragId;
  ECOSTRING chrName;
  ECOSTRING hitChrName;
  ECOSTRING seqType;
  ECOSTRING hitSeqType;
  int chrStart;
  int chrEnd;
  int hitChrStart;
  int hitChrEnd;
  int relOri;
};
#undef FUNCSTRUCTTYPE

SyntenyRegion *SyntenyRegion_new();

#define SyntenyRegion_setEnd(sr,e) (sr)->end = (e)
#define SyntenyRegion_getEnd(sr) (sr)->end

#define SyntenyRegion_setStart(sr,s) (sr)->start = (s)
#define SyntenyRegion_getStart(sr) (sr)->start

#define SyntenyRegion_setRelOri(sr,ro) (sr)->relOri = (ro)
#define SyntenyRegion_getRelOri(sr) (sr)->relOri

#define SyntenyRegion_setClusterId(sr,cid) (sr)->clusterId = (cid)
#define SyntenyRegion_getClusterId(sr) (sr)->clusterId

#define SyntenyRegion_setDNAFragId(sr,did) (sr)->dnaFragId = (did)
#define SyntenyRegion_getDNAFragId(sr) (sr)->dnaFragId

#define SyntenyRegion_setDbID(sr,dbid) Storable_setDbID(&((sr)->st),(dbid))
#define SyntenyRegion_getDbID(sr) Storable_getDbID(&((sr)->st))

#define SyntenyRegion_setAdaptor(sr,ad) Storable_setAdaptor((sr),(ad))
#define SyntenyRegion_getAdaptor(sr) Storable_getAdaptor((sr))

#define SyntenyRegion_setChrEnd(sr,e) (sr)->chrEnd = (e)
#define SyntenyRegion_getChrEnd(sr) (sr)->chrEnd

#define SyntenyRegion_setChrStart(sr,s) (sr)->chrStart = (s)
#define SyntenyRegion_getChrStart(sr) (sr)->chrStart

ECOSTRING SyntenyRegion_setChrName(SyntenyRegion *sr, char *chrName);
#define SyntenyRegion_getChrName(sr) (sr)->chrName

#define SyntenyRegion_setHitChrEnd(sr,e) (sr)->hitChrEnd = (e)
#define SyntenyRegion_getHitChrEnd(sr) (sr)->hitChrEnd

#define SyntenyRegion_setHitChrStart(sr,s) (sr)->hitChrStart = (s)
#define SyntenyRegion_getHitChrStart(sr) (sr)->hitChrStart

ECOSTRING SyntenyRegion_setHitChrName(SyntenyRegion *sr, char *hitChrName);
#define SyntenyRegion_getHitChrName(sr) (sr)->hitChrName

ECOSTRING SyntenyRegion_setHitSeqType(SyntenyRegion *sr, char *hitSeqType);

ECOSTRING SyntenyRegion_setSeqType(SyntenyRegion *sr, char *seqType);

void SyntenyRegion_free(SyntenyRegion *sr);

#ifdef __SYNTENYREGION_MAIN__
  SyntenyRegionFuncs
    syntenyRegionFuncs = {
                    SyntenyRegion_free
                   };
#else
  extern SyntenyRegionFuncs syntenyRegionFuncs;
#endif

#endif
