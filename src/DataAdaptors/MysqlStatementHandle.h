#ifndef __MYSQLSTATEMENTHANDLE_H__
#define __MYSQLSTATEMENTHANDLE_H__

#include "mysql.h"
#include "StatementHandle.h"
#include "MysqlResultRow.h"

typedef struct MysqlStatementHandleStruct MysqlStatementHandle;

#ifdef __hpux
 #include <varargs.h>
 void MysqlStatementHandle_execute(va_alist);
#else
 #include <stdarg.h>
 void MysqlStatementHandle_execute(StatementHandle *sth, ...);
#endif

StatementHandle *MysqlStatementHandle_new(DBConnection *dbc, char *query);
ResultRow *MysqlStatementHandle_fetchRow(StatementHandle *sth);
IDType MysqlStatementHandle_getInsertId(StatementHandle *sth);
void MysqlStatementHandle_finish(StatementHandle *sth);

OBJECTFUNC_TYPES(MysqlStatementHandle)

typedef struct MysqlStatementHandleFuncsStruct {
  OBJECTFUNCS_DATA(MysqlStatementHandle)
} MysqlStatementHandleFuncs;



#define MYSQLSTATEMENTHANDLE_DATA \
  STATEMENTHANDLE_DATA \
  MYSQL_RES *results; \
  MysqlResultRow *m_row;

#define FUNCSTRUCTTYPE MysqlStatementHandleFuncs
struct MysqlStatementHandleStruct {
  MYSQLSTATEMENTHANDLE_DATA
};
#undef FUNCSTRUCTTYPE

#ifdef __MYSQLSTATEMENTHANDLE_MAIN__
  MysqlStatementHandleFuncs
    mysqlStatementHandleFuncs = {
                            NULL
                           };
#else
  extern MysqlStatementHandleFuncs mysqlStatementHandleFuncs;
#endif

#endif
