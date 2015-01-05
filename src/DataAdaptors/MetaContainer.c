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

#include "MetaContainer.h"
#include "StrUtil.h"
#include "BaseAdaptor.h"
#include "MysqlUtil.h"
#include "StatementHandle.h"
#include "Species.h"

MetaContainer *MetaContainer_new(DBAdaptor *dba) {
  MetaContainer *mc;

  if ((mc = (MetaContainer *)calloc(1,sizeof(MetaContainer))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for MetaContainer\n");
    return NULL;
  }
  BaseAdaptor_init((BaseAdaptor *)mc, dba, META_CONTAINER);

  return mc;
}

Vector *MetaContainer_listValueByKey(MetaContainer *mc, char *key) {
  Vector *results = Vector_new();
  StatementHandle *sth;
  ResultRow *row;
  char qStr[1024];

  sprintf(qStr,
       "SELECT meta_value"
       " FROM meta"
       " WHERE meta_key = '%s' ORDER BY meta_id",key);

  sth = mc->prepare((BaseAdaptor *)mc,qStr,strlen(qStr)); 
  sth->execute(sth);

  while ((row = sth->fetchRow(sth))) {
    char *tmpStr;
    Vector_addElement(results, StrUtil_copyString(&tmpStr,row->getStringAt(row,0),0));
  }
  sth->finish(sth);

  return results;
}

int MetaContainer_getIntValueByKey(MetaContainer *mc, char *key, int *retVal) {
  StatementHandle *sth;
  ResultRow *row;
  char qStr[1024];
  int okFlag = 0;

  sprintf(qStr,
       "SELECT meta_value"
       " FROM meta"
       " WHERE meta_key = '%s' ORDER BY meta_id",key);

  sth = mc->prepare((BaseAdaptor *)mc,qStr,strlen(qStr)); 
  sth->execute(sth);

  if ((row = sth->fetchRow(sth))) {
    *retVal = row->getIntAt(row,0); 
    okFlag = 1;
  }
  sth->finish(sth);

  return okFlag;
}

Species *MetaContainer_getSpecies(MetaContainer *mc) {
  Vector *classification;
  StatementHandle *sth;
  ResultRow *row;
  char *commonName;
  char qStr[256];
  Species *species;

  sprintf(qStr,
             "SELECT meta_value"
             " FROM meta"
             " WHERE meta_key = 'species.common_name'");

  sth = mc->prepare((BaseAdaptor *)mc, qStr, strlen(qStr));
  sth->execute(sth);
  if ((row = sth->fetchRow(sth))) {
    commonName = row->getStringAt(row,0);
  } else {
    return NULL;
  }

  classification = MetaContainer_listValueByKey(mc, "species.classification");
  if (!classification) {
    return NULL;
  }

  species = Species_new();
  Species_setCommonName(species, commonName );
  Species_setClassification(species, classification);

  return species;

}

char *MetaContainer_getDefaultAssembly(MetaContainer *mc) {
  char *result = NULL;
  char *qStr = "SELECT meta_value from meta where meta_key = 'assembly.default'";
  StatementHandle *sth = mc->prepare((BaseAdaptor *)mc,qStr,strlen(qStr));

  if (sth) {
    sth->execute(sth);
    ResultRow *row = sth->fetchRow(sth);
    char *assStr;

    if (row) {
      assStr = row->getStringCopyAt(row,0);
      sth->finish(sth);
      result = assStr;
    } else {
      sth->finish(sth);
      result = NULL;
    }
  }

  return result;
}
