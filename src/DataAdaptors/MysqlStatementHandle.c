#include "MysqlStatementHandle.h"
#include "MysqlResultRow.h"
#include "StrUtil.h"
#include "mysql.h"
#include "EnsC.h"


StatementHandle *MysqlStatementHandle_new(DBConnection *dbc, char *query) {
  MysqlStatementHandle *sth;

  if ((sth = (MysqlStatementHandle *)calloc(1,sizeof(MysqlStatementHandle))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for sth\n");
    return NULL;
  }

  sth->objectType = CLASS_MYSQLSTATEMENTHANDLE;

  sth->execute     = MysqlStatementHandle_execute;
  sth->fetchRow    = MysqlStatementHandle_fetchRow;
  sth->finish      = MysqlStatementHandle_finish;
  sth->getInsertId = MysqlStatementHandle_getInsertId;

  sth->dbc = dbc;

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
void MysqlStatementHandle_execute(StatementHandle *sth, ...) {
  va_list args;
  char statement[EXTREMELEN];
  int qlen;
  MYSQL_RES *results;
  unsigned int field_count;
  MysqlStatementHandle *m_sth;

  Class_assertType(CLASS_MYSQLSTATEMENTHANDLE,sth->objectType);

  printf("Statement = %s\n",sth->statementFormat);

  m_sth = (MysqlStatementHandle *)sth;

  va_start(args, sth);
  qlen = vsnprintf(statement, EXTREMELEN, m_sth->statementFormat, args);
  va_end(args);

  if (qlen < 0 || qlen > EXTREMELEN) {
    fprintf(stderr, "ERROR: vsnprintf call failed during statement execution\n");
    exit(1);
  }


  if (mysql_real_query (m_sth->dbc->mysql, statement, qlen) != 0) {    /* the query failed */
    fprintf(stderr, "Could not execute query %s\n", statement);
    return;
  }

  /* the query succeeded; determine whether or not it returns data */

  results = mysql_store_result (m_sth->dbc->mysql);
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
    }
  }
}

#endif


ResultRow *MysqlStatementHandle_fetchRow(StatementHandle *sth) {
  MysqlStatementHandle *m_sth;
  MysqlResultRow *m_row;
  MYSQL_ROW mysql_row;

  Class_assertType(CLASS_MYSQLSTATEMENTHANDLE,sth->objectType);

  m_sth = (MysqlStatementHandle *)sth;

  if (m_sth->results == NULL) {
    fprintf(stderr,"ERROR: Tried to fetch a row for a StatementHandle with no results for %s\n",
            sth->statementFormat);
    exit(1);
  }

  mysql_row = mysql_fetch_row(m_sth->results);

  if (mysql_row == NULL) {
    return NULL;
  }
  m_row = MysqlResultRow_new();
  m_row->mysql_row = mysql_row;

  return (ResultRow *)m_row;
}

int64 MysqlStatementHandle_getInsertId(StatementHandle *sth) {
  MysqlStatementHandle *m_sth;
  int64 insertId;

  Class_assertType(CLASS_MYSQLSTATEMENTHANDLE,sth->objectType);

  m_sth = (MysqlStatementHandle *)sth;

  insertId = mysql_insert_id(m_sth->dbc->mysql);

  if (insertId == 0) {
    fprintf(stderr, "Warning: Insert id was 0\n");
  }

  return insertId;
}

void MysqlStatementHandle_finish(StatementHandle *sth) {
  MysqlStatementHandle *m_sth;

  Class_assertType(CLASS_MYSQLSTATEMENTHANDLE,sth->objectType);

  m_sth = (MysqlStatementHandle *)sth;
}
