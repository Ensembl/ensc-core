#ifndef __MYSQLUTIL_H__
#define __MYSQLUTIL_H__

#include <mysql.h>

int MysqlUtil_getInt(MYSQL_ROW row, int col);
long MysqlUtil_getLong(MYSQL_ROW row, int col);
double MysqlUtil_getDouble(MYSQL_ROW row, int col);
char *MysqlUtil_getString(MYSQL_ROW row, int col);

#endif
