#ifndef __MYSQLRESULTROW_H__
#define __MYSQLRESULTROW_H__

#include <mysql.h>
#include "ResultRow.h"

typedef struct MysqlResultRowStruct MysqlResultRow;

MysqlResultRow *MysqlResultRow_new();

char *    MysqlResultRow_getStringAt(ResultRow *row, int ind);
int       MysqlResultRow_getIntAt(ResultRow *row, int ind);
long      MysqlResultRow_getLongAt(ResultRow *row, int ind);
IDType     MysqlResultRow_getLongLongAt(ResultRow *row, int ind);
double    MysqlResultRow_getDoubleAt(ResultRow *row, int ind);
char *    MysqlResultRow_col(ResultRow *row, int ind);


#define MYSQLRESULTROW_DATA \
  RESULTROW_DATA \
  MYSQL_ROW mysql_row;

struct MysqlResultRowStruct {
  MYSQLRESULTROW_DATA
};

#endif
