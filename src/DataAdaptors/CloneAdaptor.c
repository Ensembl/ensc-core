#include "CloneAdaptor.h"
#include "BaseAdaptor.h"
#include "MysqlUtil.h"


CloneAdaptor *CloneAdaptor_new(DBAdaptor *dba) {
  CloneAdaptor *ca;

  if ((ca = (CloneAdaptor *)calloc(1,sizeof(CloneAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for CloneAdaptor\n");
    return NULL;
  }
  BaseAdaptor_init((BaseAdaptor *)ca, dba, CLONE_ADAPTOR);

  return ca;
}

Clone *CloneAdaptor_fetchByDbID(CloneAdaptor *ca, long dbID) {
  Clone *clone;
  char qStr[256];
  MYSQL_RES *results;
  MYSQL_ROW row;

  sprintf(qStr,
    "SELECT clone_id, logic_name,"
    "       program, program_version, program_file,"
    "       db, db_version, db_file,"
    "       module, module_version,"
    "       gff_source, gff_feature,"
    "       created, parameters"
    " FROM   clone"
    " WHERE  clone_id = %d", dbID);

  results = ca->prepare((BaseAdaptor *)ca,qStr,strlen(qStr));

  row = mysql_fetch_row(results);
  if( row == NULL ) {
    return NULL;
  }

  clone = CloneAdaptor_cloneFromRow(ca, row);

  return clone;
}

Clone *CloneAdaptor_cloneFromRow(CloneAdaptor *ca, MYSQL_ROW row) {
  Clone *clone = Clone_new();

  return clone;
}
