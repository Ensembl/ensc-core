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

#define __MYSQLSTATEMENTHANDLE_MAIN__
#include "MysqlStatementHandle.h"
#undef __MYSQLSTATEMENTHANDLE_MAIN__
#include "MysqlResultRow.h"
#include "StrUtil.h"
#include "mysql.h"
#include "EnsC.h"

#include "Error.h"
#include "Class.h"

#include "ProcUtil.h"

#include <string.h>


StatementHandle *MysqlStatementHandle_new(DBConnection *dbc, char *query) {
  MysqlStatementHandle *sth;

  if ((sth = (MysqlStatementHandle *)calloc(1,sizeof(MysqlStatementHandle))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for sth\n");
    return NULL;
  }

  sth->objectType = CLASS_MYSQLSTATEMENTHANDLE;

  sth->funcs = &mysqlStatementHandleFuncs;

  sth->execute     = MysqlStatementHandle_execute;
  sth->fetchRow    = MysqlStatementHandle_fetchRow;
  sth->numRows     = MysqlStatementHandle_numRows;
  sth->finish      = MysqlStatementHandle_finish;
  sth->getInsertId = MysqlStatementHandle_getInsertId;
  sth->addFlag     = MysqlStatementHandle_addFlag;

  sth->dbc = dbc;
 
  sth->m_row = MysqlResultRow_new();

  if ((sth->statementFormat = StrUtil_copyString(&(sth->statementFormat),
                                                 query,0)) == NULL) {
    Error_trace("MysqlStatementHandle_new", NULL);
    return NULL;
  }

  return (StatementHandle *)sth;
}


#ifdef __hpux /* Machines with vararg vsprintf */

void MysqlStatementHandle_execute(va_alist) {
  fprintf(stderr, "ERROR: vararg version of MysqlStatementHandle_execute not implemented - Nag Steve\n");
  exit(1);
}

#else /* Machines with stdarg vsprintf */

/* bits of this came from the process_query routine in the MYSQL
 * book.
 */
unsigned long long MysqlStatementHandle_execute(StatementHandle *sth, ...) {
  va_list args;
  char *statement = NULL;
  int qlen;
  MYSQL_RES *results;
  MysqlStatementHandle *m_sth;

  if ((statement = (char *)calloc(655500,sizeof(char))) == NULL) {
    fprintf(stderr,"Failed allocating statment\n");
    return 0;
  }

  Class_assertType(CLASS_MYSQLSTATEMENTHANDLE,sth->objectType);

  //printf("Statement = %s\n",sth->statementFormat);

  m_sth = (MysqlStatementHandle *)sth;

  va_start(args, sth);
  qlen = vsnprintf(statement, 655500, m_sth->statementFormat, args);
  va_end(args);

  if (qlen < 0 || qlen > 655500) {
    fprintf(stderr, "ERROR: vsnprintf call failed during statement execution (qlen = %d)\n", qlen);
    return 0;
  }

  //fprintf(stderr, "Statement after formatting = %s\n",statement);

  if (mysql_real_query (m_sth->dbc->mysql, statement, qlen) != 0) {    /* the query failed */
    fprintf(stderr, "Could not execute query %s\n\n", statement);
    fprintf(stderr, "Query length is %d\n\n", (int)strlen(statement));
    //fprintf(stderr, "Stack trace:\n");
    //ProcUtil_showBacktrace(EnsC_progName);
    fprintf(stderr, "Database %s host %s user %s pass %s port %d\n", 
            DBConnection_getDbName(m_sth->dbc),
            DBConnection_getHost(m_sth->dbc),
            DBConnection_getUser(m_sth->dbc),
            DBConnection_getPass(m_sth->dbc),
            DBConnection_getPort(m_sth->dbc));
    //if (m_sth->dbc->mysql->net && m_sth->dbc->mysql->net->last_error)
    //  fprintf(stderr, "Mysql error: %s", m_sth->dbc->mysql->net->last_error);
    
    return 0;
  }

  /* the query succeeded; determine whether or not it returns data */

  if (m_sth->flags & MYSQLFLAG_USE_RESULT) {
    results = mysql_use_result (m_sth->dbc->mysql);
  } else {
    results = mysql_store_result (m_sth->dbc->mysql);
  }
  if (results) {                    /* a result set was returned */
    /* store it for future fetches */
    m_sth->results = results;

  } else { /* no result set was returned */ 
    /*
     * does the lack of a result set mean that the query didn't
     * return one, or that it should have but an error occurred?
     */
    if (!mysql_field_count(m_sth->dbc->mysql)) {
      /*
       * query generated no result set (it was not a SELECT, SHOW,
       * DESCRIBE, etc.) - do nothing
       */
    } else {   /* an error occurred */
      fprintf (stderr, "Could not retrieve result set");
      return 0;
    }
  }

  free(statement);
  return mysql_affected_rows(m_sth->dbc->mysql);
}

#endif

void MysqlStatementHandle_addFlag(StatementHandle *sth, unsigned long flag) {
  sth->flags |= flag;
}
  

ResultRow *MysqlStatementHandle_fetchRow(StatementHandle *sth) {
  ResultRow *result = NULL;
  MysqlStatementHandle *m_sth;
  MYSQL_ROW mysql_row;

  Class_assertType(CLASS_MYSQLSTATEMENTHANDLE,sth->objectType);

  m_sth = (MysqlStatementHandle *)sth;

  if (m_sth->results == NULL) {
    fprintf(stderr,"ERROR: Tried to fetch a row for a StatementHandle with no results for %s\n",
            sth->statementFormat);
  }
  else {
    mysql_row = mysql_fetch_row(m_sth->results);

    if (mysql_row != NULL) {
      //m_row = MysqlResultRow_new();
      m_sth->m_row->mysql_row = mysql_row;
      result = (ResultRow *)(m_sth->m_row);
    }
  }

  return result;
}

unsigned long long MysqlStatementHandle_numRows(StatementHandle *sth) {
  MysqlStatementHandle *m_sth;

  Class_assertType(CLASS_MYSQLSTATEMENTHANDLE,sth->objectType);

  m_sth = (MysqlStatementHandle *)sth;

  if (m_sth->results == NULL) {
    fprintf(stderr,"ERROR: Tried to fetch number of rows for a StatementHandle with no results for %s\n",
            sth->statementFormat);
    return 0;
  }

  return mysql_num_rows(m_sth->results);
}

IDType MysqlStatementHandle_getInsertId(StatementHandle *sth) {
  MysqlStatementHandle *m_sth;
  IDType insertId;

  Class_assertType(CLASS_MYSQLSTATEMENTHANDLE,sth->objectType);

  m_sth = (MysqlStatementHandle *)sth;

  insertId = mysql_insert_id(m_sth->dbc->mysql);

  if (insertId == 0) {
    //If we do a INSERT IGNORE and insertId equals 0, it just means that the row already exists
    if (strncmp(sth->statementFormat, "INSERT IGNORE", sizeof(char)*13)) {
      fprintf(stderr, "Warning: Insert id was 0\n%s\n", sth->statementFormat);
      //ProcUtil_showBacktrace(EnsC_progName);
    }
  }

  return insertId;
}

void MysqlStatementHandle_finish(StatementHandle *sth) {
  MysqlStatementHandle *m_sth;

  Class_assertType(CLASS_MYSQLSTATEMENTHANDLE,sth->objectType);

  m_sth = (MysqlStatementHandle *)sth;

  if (m_sth->results) {
//    fprintf(stderr, "Freeing results\n");
    mysql_free_result(m_sth->results);
  }

  if (m_sth->statementFormat) free(m_sth->statementFormat);
  if (m_sth->m_row) free(m_sth->m_row);

  free(m_sth);
}
