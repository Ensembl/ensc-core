#include "CloneAdaptor.h"
#include "BaseAdaptor.h"
#include "MysqlUtil.h"

#include "StatementHandle.h"
#include "ResultRow.h"


CloneAdaptor *CloneAdaptor_new(DBAdaptor *dba) {
  CloneAdaptor *ca;

  if ((ca = (CloneAdaptor *)calloc(1,sizeof(CloneAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for CloneAdaptor\n");
    return NULL;
  }
  BaseAdaptor_init((BaseAdaptor *)ca, dba, CLONE_ADAPTOR);

  return ca;
}

Clone *CloneAdaptor_fetchByDbID(CloneAdaptor *ca, int64 dbID) {
  Clone *clone;
  char qStr[256];
  StatementHandle *sth;
  ResultRow *row;

  sprintf(qStr,
    "SELECT clone_id, logic_name,"
    "       program, program_version, program_file,"
    "       db, db_version, db_file,"
    "       module, module_version,"
    "       gff_source, gff_feature,"
    "       created, parameters"
    " FROM   clone"
    " WHERE  clone_id = "
    INT64FMTSTR, dbID);

  sth = ca->prepare((BaseAdaptor *)ca,qStr,strlen(qStr));
  sth->execute(sth);

  row = sth->fetchRow(sth);
  if( row == NULL ) {
    sth->finish(sth);
    return NULL;
  }

  clone = CloneAdaptor_cloneFromRow(ca, row);
  sth->finish(sth);

  return clone;
}

Clone *CloneAdaptor_cloneFromRow(CloneAdaptor *ca, ResultRow  *row) {
  Clone *clone = Clone_new();

  return clone;
}
