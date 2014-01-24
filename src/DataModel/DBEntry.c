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

// Not sure if this one should be char * or ECOSSTRING
char *DBEntry_setInfoText(DBEntry *dbe, char *infoText) {
  if ((dbe->infoText = (char *)malloc(strlen(infoText)+1)) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for infoText\n");
    return NULL;
  }

  strcpy(dbe->infoText,infoText);

  return dbe->infoText;
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

ECOSTRING DBEntry_setDbDisplayName(DBEntry *dbe, char *dbDisplayName) {
  EcoString_copyStr(ecoSTable, &(dbe->dbDisplayName),dbDisplayName,0); 

  if (dbe->dbDisplayName == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for dbDisplayName\n");
    return NULL;
  }

  return dbe->dbDisplayName;
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

ECOSTRING DBEntry_setVersion(DBEntry *dbe, char *version) {
  EcoString_copyStr(ecoSTable, &(dbe->version),version,0); 

  if (dbe->version == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for version\n");
    return NULL;
  }

  return dbe->version;
}

ECOSTRING DBEntry_setRelease(DBEntry *dbe, char *release) {
  EcoString_copyStr(ecoSTable, &(dbe->release),release,0); 

  if (dbe->release == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for release\n");
    return NULL;
  }

  return dbe->release;
}

ECOSTRING DBEntry_setInfoType(DBEntry *dbe, char *infoType) {
  EcoString_copyStr(ecoSTable, &(dbe->infoType),infoType,0); 

  if (dbe->infoType == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for infoType\n");
    return NULL;
  }

  return dbe->infoType;
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
  if (dbe->dbDisplayName) EcoString_freeStr(ecoSTable, dbe->dbDisplayName);
  if (dbe->status) EcoString_freeStr(ecoSTable, dbe->status);
  if (dbe->version) EcoString_freeStr(ecoSTable, dbe->version);
  if (dbe->release) EcoString_freeStr(ecoSTable, dbe->release);
  if (dbe->infoType) EcoString_freeStr(ecoSTable, dbe->infoType);

  if (dbe->displayId)   free(dbe->displayId);
  if (dbe->description) free(dbe->description);
  if (dbe->primaryId)   free(dbe->primaryId);
  if (dbe->infoText)    free(dbe->infoText);

  if (dbe->synonyms) Vector_free(dbe->synonyms);

  if (dbe->idXref) IdentityXref_free(dbe->idXref);

  free(dbe);
}

