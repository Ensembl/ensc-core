#include "DBEntry.h"

DBEntry *DBEntry_new() {
  DBEntry *dbe;

  if ((dbe = (DBEntry *)calloc(1,sizeof(DBEntry))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for dbe\n");
    return NULL;
  }

  dbe->objectType = CLASS_DBENTRY;
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
    dbe->synonyms = Set_new();
  }
  Set_addElement(dbe->synonyms, syn);
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
