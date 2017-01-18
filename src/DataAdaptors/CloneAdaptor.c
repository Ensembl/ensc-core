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

#include "CloneAdaptor.h"
#include "BaseAdaptor.h"
#include "MysqlUtil.h"

#include "StatementHandle.h"
#include "ResultRow.h"


CloneAdaptor *CloneAdaptor_new(DBAdaptor *dba) {
  CloneAdaptor *ca;

  if ((ca = (CloneAdaptor *)calloc(1,sizeof(CloneAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for CloneAdaptor\n");
    return NULL;
  }
  BaseAdaptor_init((BaseAdaptor *)ca, dba, CLONE_ADAPTOR);

  return ca;
}

Clone *CloneAdaptor_fetchByDbID(CloneAdaptor *ca, IDType dbID) {
  Clone *clone;
  char qStr[256];
  StatementHandle *sth;
  ResultRow *row;

  sprintf(qStr,
    "SELECT clone_id, logic_name,"
    "       program, program_version, program_file,"
    "       db, db_version, db_file,"
    "       module, module_version,"
    "       gff_source, gff_feature,"
    "       created, parameters"
    " FROM   clone"
    " WHERE  clone_id = "
    IDFMTSTR, dbID);

  sth = ca->prepare((BaseAdaptor *)ca,qStr,strlen(qStr));
  sth->execute(sth);

  row = sth->fetchRow(sth);
  if( row == NULL ) {
    sth->finish(sth);
    return NULL;
  }

  clone = CloneAdaptor_cloneFromRow(ca, row);
  sth->finish(sth);

  return clone;
}

Clone *CloneAdaptor_cloneFromRow(CloneAdaptor *ca, ResultRow  *row) {
  Clone *clone = Clone_new();

  return clone;
}
