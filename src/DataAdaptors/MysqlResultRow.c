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

#define __MYSQLRESULTROW_MAIN__
#include "MysqlResultRow.h"
#undef __MYSQLRESULTROW_MAIN__

#include "MysqlUtil.h"
#include "Class.h"

#include <stdlib.h>

MysqlResultRow *MysqlResultRow_new() {
  MysqlResultRow *rr;

  if ((rr = (MysqlResultRow *)calloc(1,sizeof(MysqlResultRow))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for rr\n");
    return NULL;
  }

  rr->objectType = CLASS_MYSQLRESULTROW;

  rr->funcs = &mysqlResultRowFuncs;

  rr->getStringAt     = MysqlResultRow_getStringAt;
  rr->getStringAllowNullAt = MysqlResultRow_getStringAllowNullAt;
  rr->getStringCopyAt = MysqlResultRow_getStringCopyAt;
  rr->getIntAt        = MysqlResultRow_getIntAt;
  rr->getLongAt       = MysqlResultRow_getLongAt;
  rr->getLongLongAt   = MysqlResultRow_getLongLongAt;
  rr->getDoubleAt     = MysqlResultRow_getDoubleAt;
  rr->col             = MysqlResultRow_col;

  return rr;
}

char * MysqlResultRow_getStringCopyAt(ResultRow *row, int ind) {
  MysqlResultRow *m_row;

  Class_assertType(CLASS_MYSQLRESULTROW, row->objectType);

  m_row = (MysqlResultRow *)row;

  return MysqlUtil_getStringCopy(m_row->mysql_row, ind);
}

char * MysqlResultRow_getStringAllowNullAt(ResultRow *row, int ind) {
  MysqlResultRow *m_row;

  Class_assertType(CLASS_MYSQLRESULTROW, row->objectType);

  m_row = (MysqlResultRow *)row;

  return MysqlUtil_getStringAllowNull(m_row->mysql_row, ind);
}

char * MysqlResultRow_getStringAt(ResultRow *row, int ind) {
  MysqlResultRow *m_row;

  Class_assertType(CLASS_MYSQLRESULTROW, row->objectType);

  m_row = (MysqlResultRow *)row;

  return MysqlUtil_getStringNoNull(m_row->mysql_row, ind);
}

char * MysqlResultRow_col(ResultRow *row, int ind) {
  MysqlResultRow *m_row;

  Class_assertType(CLASS_MYSQLRESULTROW, row->objectType);

  m_row = (MysqlResultRow *)row;

  return m_row->mysql_row[ind];
}

int MysqlResultRow_getIntAt(ResultRow *row, int ind) {
  MysqlResultRow *m_row;

  Class_assertType(CLASS_MYSQLRESULTROW, row->objectType);

  m_row = (MysqlResultRow *)row;

  return MysqlUtil_getInt(m_row->mysql_row, ind);
}

long MysqlResultRow_getLongAt(ResultRow *row, int ind) {
  MysqlResultRow *m_row;

  Class_assertType(CLASS_MYSQLRESULTROW, row->objectType);

  m_row = (MysqlResultRow *)row;

  return MysqlUtil_getLong(m_row->mysql_row, ind);
}

IDType MysqlResultRow_getLongLongAt(ResultRow *row, int ind) {
  MysqlResultRow *m_row;

  Class_assertType(CLASS_MYSQLRESULTROW, row->objectType);

  m_row = (MysqlResultRow *)row;

  return MysqlUtil_getLongLong(m_row->mysql_row, ind);
}

double MysqlResultRow_getDoubleAt(ResultRow *row, int ind) {
  MysqlResultRow *m_row;

  Class_assertType(CLASS_MYSQLRESULTROW, row->objectType);

  m_row = (MysqlResultRow *)row;

  return MysqlUtil_getDouble(m_row->mysql_row, ind);
}
