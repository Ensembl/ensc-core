#ifndef __STATEMENTHANDLE_H__
#define __STATEMENTHANDLE_H__

#include "Object.h"
#include "ResultRow.h"
#include "DBConnection.h"

typedef struct StatementHandleStruct StatementHandle;

typedef void (*StatementHandle_ExecuteFunc)(StatementHandle *sth, ...);
typedef ResultRow *(*StatementHandle_FetchRowFunc)(StatementHandle *sth);
typedef void (*StatementHandle_FinishFunc)(StatementHandle *sth);

#define STATEMENTHANDLE_DATA \
  OBJECT_DATA \
  char *statementFormat; \
  char *currentStatement; \
  DBConnection *dbc; \
  StatementHandle_ExecuteFunc execute; \
  StatementHandle_FetchRowFunc fetchRow; \
  StatementHandle_FinishFunc finish;
  
struct StatementHandleStruct {
  STATEMENTHANDLE_DATA
};

#endif
