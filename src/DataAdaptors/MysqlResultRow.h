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

#ifndef __MYSQLRESULTROW_H__
#define __MYSQLRESULTROW_H__

#include <mysql.h>
#include "ResultRow.h"

typedef struct MysqlResultRowStruct MysqlResultRow;

MysqlResultRow *MysqlResultRow_new();

char *    MysqlResultRow_getStringAt(ResultRow *row, int ind);
char *    MysqlResultRow_getStringCopyAt(ResultRow *row, int ind);
char *    MysqlResultRow_getStringAllowNullAt(ResultRow *row, int ind);
int       MysqlResultRow_getIntAt(ResultRow *row, int ind);
long      MysqlResultRow_getLongAt(ResultRow *row, int ind);
IDType    MysqlResultRow_getLongLongAt(ResultRow *row, int ind);
double    MysqlResultRow_getDoubleAt(ResultRow *row, int ind);
char *    MysqlResultRow_col(ResultRow *row, int ind);

OBJECTFUNC_TYPES(MysqlResultRow)

typedef struct MysqlResultRowFuncsStruct {
  OBJECTFUNCS_DATA(MysqlResultRow)
} MysqlResultRowFuncs;


#define MYSQLRESULTROW_DATA \
  RESULTROW_DATA \
  MYSQL_ROW mysql_row;

#define FUNCSTRUCTTYPE MysqlResultRowFuncs
struct MysqlResultRowStruct {
  MYSQLRESULTROW_DATA
};
#undef FUNCSTRUCTTYPE

#ifdef __MYSQLRESULTROW_MAIN__
  MysqlResultRowFuncs
    mysqlResultRowFuncs = {
                      NULL, // free
                      NULL, // shallowCopy
                      NULL  // deepCopy
                     };
#else
  extern MysqlResultRowFuncs mysqlResultRowFuncs;
#endif


#endif
