#ifndef __DBENTRY_H__
#define __DBENTRY_H__

#include "DataModelTypes.h"
#include "Storable.h"
#include "Vector.h"
#include "EcoString.h"
#include "IdentityXref.h"
#include "Object.h"

struct DBEntryStruct {
  OBJECT_DATA
  Storable  st;
  char     *primaryId;
  ECOSTRING dbName;
  ECOSTRING status;
  int       version;
  char     *displayId;
  int       release;
  Vector   *synonyms;
  char     *description;
  IdentityXref *idXref;
};

DBEntry *DBEntry_new(void);

#define DBEntry_setDbID(d,dbID) Storable_setDbID(&((d)->st),dbID)
#define DBEntry_getDbID(d) Storable_getDbID(&((d)->st))

#define DBEntry_setAdaptor(d,ad) Storable_setAdaptor(&((d)->st),ad)
#define DBEntry_getAdaptor(d) Storable_getAdaptor(&((d)->st))

#define DBEntry_setVersion(d,ver) (d)->version = (ver)
#define DBEntry_getVersion(d) (d)->version

#define DBEntry_setRelease(d,rel) (d)->release = (rel)
#define DBEntry_getRelease(d) (d)->release

#define DBEntry_setIdentityXref(d,idx) (d)->idXref = (idx)
#define DBEntry_getIdentityXref(d) (d)->idXref

ECOSTRING DBEntry_setDbName(DBEntry *dbe, char *name);
#define DBEntry_getDbName(d) (d)->dbName

char *DBEntry_setPrimaryId(DBEntry *dbe, char *id);
#define DBEntry_getPrimaryId(d) (d)->primaryId

char *DBEntry_setDisplayId(DBEntry *dbe, char *id);
#define DBEntry_getDisplayId(d) (d)->displayId

char *DBEntry_setDescription(DBEntry *dbe, char *desc);
#define DBEntry_getDescription(d) (d)->description


int DBEntry_addSynonym(DBEntry *dbe, char *syn);
#define DBEntry_getAllSynonyms(d) (d)->synonyms

ECOSTRING DBEntry_setStatus(DBEntry *dbe, char *status);


#endif
