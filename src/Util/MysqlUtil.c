#include "MysqlUtil.h"
#include "StrUtil.h"

#include <stdlib.h>

int MysqlUtil_getInt(MYSQL_ROW row, int col) {
  if (row[col] == NULL) {
    return 0;
  } else {
    return atol(row[col]);
  }
}

long MysqlUtil_getLong(MYSQL_ROW row, int col) {
  if (row[col] == NULL) {
    return 0;
  } else {
    return atol(row[col]);
  }
}

IDType MysqlUtil_getLongLong(MYSQL_ROW row, int col) {
  long long val;

  if (row[col] == NULL) {
    val = 0;
  } else {
#ifdef __osf__
    if (sscanf(row[col],"%Ld",&val) == 0) {
#else
    if (sscanf(row[col],"%qd",&val) == 0) {
#endif
      fprintf(stderr,"Error: Failed to decode a long long from %s\n",row[col]);
    }
  }

  return val;
}

double MysqlUtil_getDouble(MYSQL_ROW row, int col) {
  return atof(row[col]);
}

// Doesn't make a copy of string
char *MysqlUtil_getString(MYSQL_ROW row, int col) {
  char *copy;
  if (row[col] == NULL) {
    return "";
  } else {
    return row[col];
  }
}

// Makes a copy of string
char *MysqlUtil_getStringCopy(MYSQL_ROW row, int col) {
  char *copy;
  if (row[col] == NULL) {
    if ((copy = StrUtil_copyString(&copy,"",0)) == NULL) {
      fprintf(stderr,"ERROR: Failed copying mysql col\n");
      return NULL;
    }
  } else {
    if ((copy = StrUtil_copyString(&copy,row[col],0)) == NULL) {
      fprintf(stderr,"ERROR: Failed copying mysql col\n");
      return NULL;
    }
  }
  return copy;
}
