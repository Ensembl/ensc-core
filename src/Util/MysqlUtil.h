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
