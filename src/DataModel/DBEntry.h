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

#ifndef __DBENTRY_H__
#define __DBENTRY_H__

#include "DataModelTypes.h"
#include "Storable.h"
#include "Vector.h"
#include "EcoString.h"
#include "IdentityXref.h"
#include "Object.h"

OBJECTFUNC_TYPES(DBEntry)

typedef struct DBEntryFuncsStruct {
  OBJECTFUNCS_DATA(DBEntry)
} DBEntryFuncs;

#define FUNCSTRUCTTYPE DBEntryFuncs
struct DBEntryStruct {
  OBJECT_DATA
  Storable  st;
  char     *primaryId;
  ECOSTRING dbName;
  ECOSTRING dbDisplayName;
  ECOSTRING status;
  ECOSTRING version;
  char     *displayId;
  ECOSTRING release;
  ECOSTRING infoType;
  char     *infoText;
  Vector   *synonyms;
  char     *description;
  IdentityXref *idXref;
};
#undef FUNCSTRUCTTYPE


DBEntry *DBEntry_new(void);

#define DBEntry_isStored(d, db) Storable_isStored(&((d)->st), (db))

#define DBEntry_setDbID(d,dbID) Storable_setDbID(&((d)->st),dbID)
#define DBEntry_getDbID(d) Storable_getDbID(&((d)->st))

#define DBEntry_setAdaptor(d,ad) Storable_setAdaptor(&((d)->st),ad)
#define DBEntry_getAdaptor(d) Storable_getAdaptor(&((d)->st))

//#define DBEntry_setVersion(d,ver) (d)->version = (ver)
ECOSTRING DBEntry_setVersion(DBEntry *dbe, char *version);
#define DBEntry_getVersion(d) (d)->version

//#define DBEntry_setRelease(d,rel) (d)->release = (rel)
ECOSTRING DBEntry_setRelease(DBEntry *dbe, char *release);
#define DBEntry_getRelease(d) (d)->release

#define DBEntry_setIdentityXref(d,idx) (d)->idXref = (idx)
#define DBEntry_getIdentityXref(d) (d)->idXref

ECOSTRING DBEntry_setDbName(DBEntry *dbe, char *name);
#define DBEntry_getDbName(d) (d)->dbName

ECOSTRING DBEntry_setDbDisplayName(DBEntry *dbe, char *name);
#define DBEntry_getDbDisplayName(d) (d)->dbDisplayName

char *DBEntry_setPrimaryId(DBEntry *dbe, char *id);
#define DBEntry_getPrimaryId(d) (d)->primaryId

char *DBEntry_setDisplayId(DBEntry *dbe, char *id);
#define DBEntry_getDisplayId(d) (d)->displayId

char *DBEntry_setDescription(DBEntry *dbe, char *desc);
#define DBEntry_getDescription(d) (d)->description

int DBEntry_addSynonym(DBEntry *dbe, char *syn);
#define DBEntry_getAllSynonyms(d) (d)->synonyms

ECOSTRING DBEntry_setStatus(DBEntry *dbe, char *status);
#define DBEntry_getStatus(d) (d)->status

ECOSTRING DBEntry_setInfoType(DBEntry *dbe, char *infoType);
#define DBEntry_getInfoType(d) (d)->infoType

char *DBEntry_setInfoText(DBEntry *dbe, char *infoText);
#define DBEntry_getInfoText(d) (d)->infoText


void DBEntry_free(DBEntry *dbe);


#ifdef __DBENTRY_MAIN__
  DBEntryFuncs
    dBEntryFuncs = {
                    DBEntry_free,
                    NULL, // shallowCopy
                    NULL  // deepCopy
                   };
#else
  extern DBEntryFuncs dBEntryFuncs;
#endif


#endif
