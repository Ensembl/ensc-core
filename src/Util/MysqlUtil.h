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

#ifndef __MYSQLUTIL_H__
#define __MYSQLUTIL_H__

#include <mysql.h>
#include "EnsC.h"

int     MysqlUtil_getInt(MYSQL_ROW row, int col);
long    MysqlUtil_getLong(MYSQL_ROW row, int col);
IDType  MysqlUtil_getLongLong(MYSQL_ROW row, int col);
double  MysqlUtil_getDouble(MYSQL_ROW row, int col);
char *  MysqlUtil_getStringAllowNull(MYSQL_ROW row, int col);
char *  MysqlUtil_getStringNoNull(MYSQL_ROW row, int col);
char *  MysqlUtil_getStringCopy(MYSQL_ROW row, int col);

#endif
