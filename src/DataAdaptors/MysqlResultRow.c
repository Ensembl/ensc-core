#include "MysqlResultRow.h"

#include "MysqlUtil.h"
#include "Class.h"

MysqlResultRow *MysqlResultRow_new() {
  MysqlResultRow *rr;

  if ((rr = (MysqlResultRow *)calloc(1,sizeof(MysqlResultRow))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for rr\n");
    return NULL;
  }

  rr->objectType = CLASS_MYSQLRESULTROW;

  rr->getStringAt   = MysqlResultRow_getStringAt;
  rr->getIntAt      = MysqlResultRow_getIntAt;
  rr->getLongAt     = MysqlResultRow_getLongAt;
  rr->getLongLongAt = MysqlResultRow_getLongLongAt;
  rr->getDoubleAt   = MysqlResultRow_getDoubleAt;
  rr->col           = MysqlResultRow_col;

  return rr;
}

char * MysqlResultRow_getStringAt(ResultRow *row, int ind) {
  MysqlResultRow *m_row;

  Class_assertType(CLASS_MYSQLRESULTROW, row->objectType);

  m_row = (MysqlResultRow *)row;

  return MysqlUtil_getString(m_row->mysql_row, ind);
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
