#include "SliceAdaptor.h"
#include "BaseAdaptor.h"
#include "MysqlUtil.h"


SliceAdaptor *SliceAdaptor_new(DBAdaptor *dba) {
  SliceAdaptor *sa;

  if ((sa = (SliceAdaptor *)calloc(1,sizeof(SliceAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for SliceAdaptor\n");
    return NULL;
  }
  BaseAdaptor_init((BaseAdaptor *)sa, dba, SLICE_ADAPTOR);

  return sa;
}

Slice *SliceAdaptor_fetchByDbID(SliceAdaptor *sa, long dbID) {
  Slice *slice;
  char qStr[256];
  MYSQL_RES *results;
  MYSQL_ROW row;

  sprintf(qStr,
    "SELECT slice_id, logic_name,"
    "       program, program_version, program_file,"
    "       db, db_version, db_file,"
    "       module, module_version,"
    "       gff_source, gff_feature,"
    "       created, parameters"
    " FROM   slice"
    " WHERE  slice_id = %d", dbID);

  results = sa->prepare((BaseAdaptor *)sa,qStr,strlen(qStr));

  row = mysql_fetch_row(results);
  if( row == NULL ) {
    return NULL;
  }

  slice = SliceAdaptor_sliceFromRow(sa, row);

  return slice;
}

Slice *SliceAdaptor_sliceFromRow(SliceAdaptor *sa, MYSQL_ROW row) {
  Slice *slice = Slice_new();

  return slice;
}
