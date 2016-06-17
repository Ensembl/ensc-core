/*
 * Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

#include "RepeatConsensusAdaptor.h"
#include "BaseAdaptor.h"

#include "StatementHandle.h"
#include "ResultRow.h"


RepeatConsensusAdaptor *RepeatConsensusAdaptor_new(DBAdaptor *dba) {
  RepeatConsensusAdaptor *rca;

  if ((rca = (RepeatConsensusAdaptor *)calloc(1,sizeof(RepeatConsensusAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for RepeatConsensusAdaptor\n");
    return NULL;
  }
  BaseAdaptor_init((BaseAdaptor *)rca, dba, REPEATCONSENSUS_ADAPTOR);

  return rca;
}

RepeatConsensus *RepeatConsensusAdaptor_fetchByDbID(RepeatConsensusAdaptor *rca, IDType dbID) {
  char constraintStr[256];
  Vector *rcVector;
  RepeatConsensus *rc;
  
  sprintf(constraintStr,"repeat_consensus_id = " IDFMTSTR " limit 1", dbID);
  rcVector = RepeatConsensusAdaptor_genericFetch(rca, constraintStr); 

  rc = Vector_getElementAt(rcVector,0);

  Vector_free(rcVector);

  return rc;   
}

RepeatConsensus *RepeatConsensusAdaptor_fetchByName(RepeatConsensusAdaptor *rca, char *name) {
  char constraintStr[256];
  Vector *rcVector;
  RepeatConsensus *rc;
  
  sprintf(constraintStr,"repeat_name = \'%s\' limit 1", name);
  rcVector = RepeatConsensusAdaptor_genericFetch(rca, constraintStr); 

  rc = Vector_getElementAt(rcVector,0);

  Vector_free(rcVector);

  return rc;   
}

RepeatConsensus *RepeatConsensusAdaptor_fetchByNameAndClass(RepeatConsensusAdaptor *rca, char *name, char *class) {
  char constraintStr[256];
  Vector *rcVector;
  RepeatConsensus *rc;
  
  sprintf(constraintStr,"repeat_name = \'%s\' AND repeat_class = \'%s\' limit 1", name,class);
  rcVector = RepeatConsensusAdaptor_genericFetch(rca, constraintStr); 

  rc = Vector_getElementAt(rcVector,0);

  Vector_free(rcVector);

  return rc;   
}


Vector *RepeatConsensusAdaptor_fetchByClassAndSeq(RepeatConsensusAdaptor *rca, char *class, char *seq) {
  Vector *result = NULL;
  char *constraintStr = NULL;

  if ((constraintStr = (char *)calloc(655500,sizeof(char))) == NULL) {
    fprintf(stderr,"Failed allocating constraintStr\n");
    return result;
  }
  
  sprintf(constraintStr,"repeat_class = \'%s\' AND repeat_consensus = \'%s\'", class, seq);
  result = RepeatConsensusAdaptor_genericFetch(rca, constraintStr); 
  
  free(constraintStr);
  return result;
}

Vector *RepeatConsensusAdaptor_genericFetch(RepeatConsensusAdaptor *rca, char *whereClause) {
  StatementHandle *sth;
  ResultRow *row;
  char *qStr = NULL;
  Vector *consensi;
    
  if ((qStr = (char *)calloc(655500,sizeof(char))) == NULL) {
    fprintf(stderr,"Failed allocating qStr\n");
    return NULL;
  }

  sprintf(qStr,"SELECT repeat_consensus_id, repeat_name,"
               "       repeat_class, LENGTH(repeat_consensus)"
               " FROM repeat_consensus"
               " WHERE %s", whereClause);

  sth = rca->prepare((BaseAdaptor *)rca,qStr,strlen(qStr));
  sth->execute(sth);

  consensi = Vector_new();
  
  while ((row = sth->fetchRow(sth))) {
    RepeatConsensus *rc = RepeatConsensus_new();

    RepeatConsensus_setDbID(rc, row->getLongLongAt(row,0));
    RepeatConsensus_setName(rc, row->getStringAt(row,1));
    RepeatConsensus_setRepeatClass(rc, row->getStringAt(row,2));
    // Do I really need to do this???? RepeatConsensus_setLength(rc, row->getIntAt(row,3));
    RepeatConsensus_setAdaptor(rc, (BaseAdaptor *)rca);

    Vector_addElement(consensi, rc);
  }

  free(qStr);
  return consensi;
}

int RepeatConsensusAdaptor_store(RepeatConsensusAdaptor *rca, Vector *consensi) {
  StatementHandle *sth;
  char qStr[1024];
  int i;

  sprintf(qStr,
    "INSERT into repeat_consensus( repeat_consensus_id"
          ", repeat_name"
          ", repeat_class"
          ", repeat_type"
          ", repeat_consensus )"
      "VALUES (NULL, '%%s','%%s','%%s','%%s')");

  sth = rca->prepare((BaseAdaptor *)rca,qStr,strlen(qStr));
    
  for (i=0; i<Vector_getNumElement(consensi); i++) {
    RepeatConsensus *rc = Vector_getElementAt(consensi,i);
    IDType dbID; 
    char *name;
    char *class;
    char *type;
    char *seq;

    if (!(name  = RepeatConsensus_getName(rc))) {
      fprintf(stderr,"Error: Name not set for repeat consensus\n");
      exit(1);
    }
    if (!(class  = RepeatConsensus_getRepeatClass(rc))) {
      fprintf(stderr,"Error: Class not set for repeat consensus\n");
      exit(1);
    }
    if (!(type  = RepeatConsensus_getRepeatType(rc))) {
      type = "";
    }
    if (!(seq  = RepeatConsensus_getConsensus(rc))) {
      fprintf(stderr,"Error: Consensus seq not set for repeat consensus\n");
      exit(1);
    }
  
    sth->execute(sth, name, class, type, seq);
    
    dbID = sth->getInsertId(sth);
    
    RepeatConsensus_setDbID(rc, dbID);
    RepeatConsensus_setAdaptor(rc, (BaseAdaptor *)rca);
  }
  sth->finish(sth);
  return 1;
}

int RepeatConsensusAdaptor_free(RepeatConsensusAdaptor *rc) {
  fprintf(stderr,"RepeatConsensusAdaptor_free not implemented\n");
  return 0;
}
