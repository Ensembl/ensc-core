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
  double val;
  if (row[col] == NULL) {
    val = 0.0;
  } else {
    val = atof(row[col]);
  }
  return val;
}

// Doesn't make a copy of string
char *MysqlUtil_getStringNoNull(MYSQL_ROW row, int col) {
  if (row[col] == NULL) {
    return "";
  } else {
    return row[col];
  }
}

char *MysqlUtil_getStringAllowNull(MYSQL_ROW row, int col) {
  if (row[col] == NULL) {
    return NULL;
  } else {
    return row[col];
  }
}

// Makes a copy of string
char *MysqlUtil_getStringCopy(MYSQL_ROW row, int col) {
  char *copy;
  if (row[col] == NULL) {
    if ((StrUtil_copyString(&copy,"",0)) == NULL) {
      fprintf(stderr,"ERROR: Failed copying mysql col\n");
      return NULL;
    }
  } else {
    if ((StrUtil_copyString(&copy,row[col],0)) == NULL) {
      fprintf(stderr,"ERROR: Failed copying mysql col\n");
      return NULL;
    }
  }
// Sanity check
  if (copy == row[col]) {
    fprintf(stderr,"ERROR: copy == row[col] in copying mysql col\n");
    exit(1);
  }
    
  return copy;
}
