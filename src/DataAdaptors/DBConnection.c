#include "DBConnection.h"
#include "StrUtil.h"

DBConnection *DBConnection_new(char *host, char *user, char *pass, 
                               char *dbname, unsigned int port) {
  DBConnection *dbc;
  MYSQL *mysql;

  if ((dbc = (DBConnection *)calloc(1,sizeof(DBConnection))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating DBConnection\n");
    return NULL;
  }

  mysql = mysql_init(NULL);

  if ((mysql = mysql_real_connect(mysql,host, user, pass, dbname, port, NULL, 0)) == NULL) {
    fprintf(stderr,"ERROR: Failed connecting to database %s (host %s user %s pass %s port %d)\n",
            dbname,host,user,pass,port);
    return NULL;
  }

  dbc->host    = StrUtil_CopyString(host);
  dbc->user    = StrUtil_CopyString(user);
  if (pass) dbc->pass    = StrUtil_CopyString(pass);
  dbc->dbname  = StrUtil_CopyString(dbname);
  dbc->mysql   = mysql;
  dbc->prepare = DBConnection_prepare;
 
  return dbc;
}

MYSQL_RES *DBConnection_prepare(DBConnection *dbc, char *queryStr, int queryLen) {
  /*fprintf(stderr,"Query = %s\n",queryStr);*/
  mysql_real_query(dbc->mysql, queryStr, queryLen);
  return mysql_store_result(dbc->mysql);
}

BaseAdaptor *DBConnection_getAdaptor(DBConnection *dbc, int type) {
  int i;
  BaseAdaptor *ad = NULL;

  for (i=0; i<dbc->nAdaptor && !ad; i++) {
    if (dbc->adaptors[i]->adaptorType == type) {
      ad = dbc->adaptors[i];
    }
  }
  if (!ad) {
    fprintf(stderr,"ERROR: No adaptor of type %s\n",Adaptor_TypeStrings[type]);
  }
  return ad;
}

int DBConnection_addAdaptor(DBConnection *dbc, BaseAdaptor *ba) {
  if (!dbc->nAdaptor || !(dbc->nAdaptor%10)) {
    if ((dbc->adaptors = (BaseAdaptor **)realloc(dbc->adaptors,(dbc->nAdaptor+10)*sizeof(BaseAdaptor *))) == NULL) {
      fprintf(stderr,"ERROR: Failed allocating adaptor array\n");
      return 0;
    }
  }

  dbc->adaptors[dbc->nAdaptor++] = ba;

  return 1;
}
