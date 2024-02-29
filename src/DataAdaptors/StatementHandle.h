/*
 * Copyright [1999-2024] EMBL-European Bioinformatics Institute
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

#ifndef __STATEMENTHANDLE_H__
#define __STATEMENTHANDLE_H__

#include "Object.h"
#include "ResultRow.h"
#include "DBConnection.h"

typedef unsigned long long (*StatementHandle_ExecuteFunc)(StatementHandle *sth, ...);
typedef ResultRow *(*StatementHandle_FetchRowFunc)(StatementHandle *sth);
typedef void (*StatementHandle_FinishFunc)(StatementHandle *sth);
typedef unsigned long long (*StatementHandle_NumRowsFunc)(StatementHandle *sth);
typedef IDType (*StatementHandle_GetInsertIdFunc)(StatementHandle *sth);
typedef void (*StatementHandle_addFlagFunc)(StatementHandle *sth, unsigned long flag);


OBJECTFUNC_TYPES(StatementHandle)

typedef struct StatementHandleFuncsStruct {
  OBJECTFUNCS_DATA(StatementHandle)
} StatementHandleFuncs;


#define STATEMENTHANDLE_DATA \
  OBJECT_DATA \
  char *statementFormat; \
  char *currentStatement; \
  DBConnection *dbc; \
  StatementHandle_ExecuteFunc execute; \
  StatementHandle_FetchRowFunc fetchRow; \
  StatementHandle_NumRowsFunc numRows; \
  StatementHandle_FinishFunc finish; \
  StatementHandle_GetInsertIdFunc getInsertId; \
  StatementHandle_addFlagFunc addFlag; \
  unsigned long flags;
  
#define FUNCSTRUCTTYPE StatementHandleFuncs
struct StatementHandleStruct {
  STATEMENTHANDLE_DATA
};
#undef FUNCSTRUCTTYPE

#endif
