#include "MysqlStatementHandle.h"
#include "mysql.h"
#include "EnsC.h"


MysqlStatementHandle *MysqlStatementHandle_new() {
  MysqlStatementHandle *sth;

  if ((sth = (MysqlStatementHandle *)calloc(1,sizeof(MysqlStatementHandle))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for sth\n");
    return NULL;
  }

  sth->objectType = CLASS_MYSQLSTATEMENTHANDLE;

  sth->execute  = MysqlStatementHandle_execute;
  sth->fetchRow = MysqlStatementHandle_fetchRow;
  sth->finish   = MysqlStatementHandle_finish;

  return sth;
}

#ifdef __hpux /* Machines with vararg vsprintf */

void MysqlStatementHandle_execute(va_alist) {
  fprintf(stderr, "ERROR: vararg version of MysqlStatementHandle_execute not implemented - Nag Steve\n");
  exit(1);
}

#else /* Machines with stdarg vsprintf */

void MysqlStatementHandle_execute(StatementHandle *sth, ...) {
  va_list args;
  char statement[EXTREMELEN];
  int retval;
  MYSQL_RES *results;

  va_start(args, sth);
  retval = vsnprintf(statement, EXTREMELEN, sth->statementFormat, args);
  va_end(args);

  if (retval < 0 || retval > EXTREMELEN) {
    fprintf(stderr, "ERROR: vsnprintf call failed during statement execution\n");
    exit(1);
  }
}

#endif


ResultRow *MysqlStatementHandle_fetchRow(StatementHandle *sth) {
}

void MysqlStatementHandle_finish(StatementHandle *sth) {
}
