#include "MysqlUtil.h"
#include "StrUtil.h"

int MysqlUtil_getInt(MYSQL_ROW row, int col) {
  return atol(row[col]);
}

long MysqlUtil_getLong(MYSQL_ROW row, int col) {
  return atol(row[col]);
}

double MysqlUtil_getDouble(MYSQL_ROW row, int col) {
  return atof(row[col]);
}

char *MysqlUtil_getString(MYSQL_ROW row, int col) {
  char *copy;
  if ((copy = StrUtil_copyString(&copy,row[col],0)) == NULL) {
    fprintf(stderr,"ERROR: Failed copying mysql col\n");
    return NULL;
  }
  return copy;
}
