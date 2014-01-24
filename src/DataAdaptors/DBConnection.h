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

#ifndef __DBCONNECTION_H__
#define __DBCONNECTION_H__

//#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "StatementHandle.h"
#include "EcoString.h"

#include <mysql.h>


typedef StatementHandle *(*DBConnection_PrepareFunc)(DBConnection *dbc, char *queryStr, int queryLen);

struct DBConnectionStruct {
  ECOSTRING host;
  ECOSTRING user;
  ECOSTRING pass;
  unsigned int   port;
  ECOSTRING dbName;
  MYSQL *mysql;
  DBConnection_PrepareFunc prepare;
  BaseAdaptor **adaptors;
  int nAdaptor;
};

DBConnection    *DBConnection_new(char *host, char *user, char *pass, char *dbname, unsigned int port);
BaseAdaptor     *DBConnection_getAdaptor(DBConnection *dbc, int type);
StatementHandle *DBConnection_prepare(DBConnection *dbc, char *queryStr, int queryLen);
int DBConnection_addAdaptor(DBConnection *dbc, BaseAdaptor *ba);
char *DBConnection_getDriverName(DBConnection *dbc);
void DBConnection_fromDateToSeconds(DBConnection *dbc, char *column, char *wrappedColumn);


#define DBConnection_getHost(dbc) (dbc)->host
#define DBConnection_getPort(dbc) (dbc)->port
#define DBConnection_getDbName(dbc) (dbc)->dbName
#define DBConnection_getPass(dbc) (dbc)->pass
#define DBConnection_getUser(dbc) (dbc)->user


#endif
