#include "DBEntry.h"

DBEntry *DBEntry_new() {
  DBEntry *dbe;

  if ((dbe = (DBEntry *)calloc(1,sizeof(DBEntry))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for dbe\n");
    return NULL;
  }

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

ECOSTRING DBEntry_setDbName(DBEntry *dbe, char *dbName) {
  EcoString_copyStr(ecoSTable, &(dbe->dbName),dbName,0); 

  if (dbe->dbName == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for dbName\n");
    return NULL;
  }

  return dbe->dbName;
}
