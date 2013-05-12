#define __DBENTRY_MAIN__
#include "DBEntry.h"
#undef __DBENTRY_MAIN__

#include <string.h>

DBEntry *DBEntry_new() {
  DBEntry *dbe;

  if ((dbe = (DBEntry *)calloc(1,sizeof(DBEntry))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for dbe\n");
    return NULL;
  }

  dbe->objectType = CLASS_DBENTRY;

  dbe->funcs = &dBEntryFuncs;

  Object_incRefCount(dbe);
  return dbe;
}

char *DBEntry_setDisplayId(DBEntry *dbe, char *displayId) {
  if ((dbe->displayId = (char *)malloc(strlen(displayId)+1)) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for displayId\n");
    return NULL;
  }

  strcpy(dbe->displayId,displayId);

  return dbe->displayId;
}

char *DBEntry_setPrimaryId(DBEntry *dbe, char *primaryId) {
  if ((dbe->primaryId = (char *)malloc(strlen(primaryId)+1)) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for primaryId\n");
    return NULL;
  }

  strcpy(dbe->primaryId,primaryId);

  return dbe->primaryId;
}

char *DBEntry_setDescription(DBEntry *dbe, char *description) {
  if ((dbe->description = (char *)malloc(strlen(description)+1)) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for description\n");
    return NULL;
  }

  strcpy(dbe->description,description);

  return dbe->description;
}

ECOSTRING DBEntry_setDbName(DBEntry *dbe, char *dbName) {
  EcoString_copyStr(ecoSTable, &(dbe->dbName),dbName,0); 

  if (dbe->dbName == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for dbName\n");
    return NULL;
  }

  return dbe->dbName;
}

int DBEntry_addSynonym(DBEntry *dbe, char *syn) {
  if (!dbe->synonyms) {
    dbe->synonyms = Vector_new();
  }
  Vector_addElement(dbe->synonyms, syn);
  return 1;
}

ECOSTRING DBEntry_setStatus(DBEntry *dbe, char *status) {
  EcoString_copyStr(ecoSTable, &(dbe->status),status,0); 

  if (dbe->status == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for status\n");
    return NULL;
  }

  return dbe->status;
}

void DBEntry_free(DBEntry *dbe) {
  Object_decRefCount(dbe);

  if (Object_getRefCount(dbe) > 0) {
    return;
  } else if (Object_getRefCount(dbe) < 0) {
    fprintf(stderr,"Error: Negative reference count for DBEntry\n"
                   "       Freeing it anyway\n");
  }

  if (dbe->dbName) EcoString_freeStr(ecoSTable, dbe->dbName);
  if (dbe->status) EcoString_freeStr(ecoSTable, dbe->status);

  if (dbe->displayId)   free(dbe->displayId);
  if (dbe->description) free(dbe->description);
  if (dbe->primaryId)   free(dbe->primaryId);

  if (dbe->synonyms) Vector_free(dbe->synonyms);

  if (dbe->idXref) IdentityXref_free(dbe->idXref);

  free(dbe);
}

