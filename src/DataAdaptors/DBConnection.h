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


#endif
