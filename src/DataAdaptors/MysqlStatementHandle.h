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

#ifndef __MYSQLSTATEMENTHANDLE_H__
#define __MYSQLSTATEMENTHANDLE_H__

#include "mysql.h"
#include "StatementHandle.h"
#include "MysqlResultRow.h"

// flag 
#define MYSQLFLAG_USE_RESULT 2

typedef struct MysqlStatementHandleStruct MysqlStatementHandle;

#ifdef __hpux
 #include <varargs.h>
 unsigned long long MysqlStatementHandle_execute(va_alist);
#else
 #include <stdarg.h>
 unsigned long long MysqlStatementHandle_execute(StatementHandle *sth, ...);
#endif

StatementHandle *MysqlStatementHandle_new(DBConnection *dbc, char *query);
ResultRow *MysqlStatementHandle_fetchRow(StatementHandle *sth);
IDType MysqlStatementHandle_getInsertId(StatementHandle *sth);
unsigned long long MysqlStatementHandle_numRows(StatementHandle *sth);
void MysqlStatementHandle_finish(StatementHandle *sth);
void MysqlStatementHandle_addFlag(StatementHandle *sth, unsigned long flag);

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
                            NULL,
                            NULL, // shallowCopy
                            NULL  // deepCopy
                           };
#else
  extern MysqlStatementHandleFuncs mysqlStatementHandleFuncs;
#endif

#endif
