#include "AnalysisAdaptor.h"
#include "BaseAdaptor.h"
#include "MysqlUtil.h"


AnalysisAdaptor *AnalysisAdaptor_new(DBAdaptor *dba) {
  AnalysisAdaptor *aa;

  if ((aa = (AnalysisAdaptor *)calloc(1,sizeof(AnalysisAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for AnalysisAdaptor\n");
    return NULL;
  }
  BaseAdaptor_init((BaseAdaptor *)aa, dba, ANALYSIS_ADAPTOR);

  aa->analCache = IDHash_new(IDHASH_SMALL);
  aa->logicNameCache = StringHash_new(STRINGHASH_SMALL);
  
  AnalysisAdaptor_fetchAll(aa);

  return aa;
}

Analysis **AnalysisAdaptor_fetchAll(AnalysisAdaptor *aa) {
  char qStr[256];
  MYSQL_RES *results;
  MYSQL_ROW row;

  sprintf(qStr,
    "SELECT analysis_id, logic_name,"
    "       program, program_version, program_file,"
    "       db, db_version, db_file,"
    "       module, module_version,"
    "       gff_source, gff_feature,"
    "       created, parameters"
    " FROM   analysis");

  results = aa->prepare((BaseAdaptor *)aa,qStr,strlen(qStr));

  while(row = mysql_fetch_row(results)) {
    Analysis *anal = AnalysisAdaptor_analysisFromRow(aa, row);

    IDHash_add(aa->analCache, Analysis_getDbID(anal), anal);
    StringHash_add(aa->logicNameCache, Analysis_getLogicName(anal), anal);
  }

  return (Analysis **)IDHash_getValues(aa->analCache);
}


Analysis *AnalysisAdaptor_fetchByDbID(AnalysisAdaptor *aa, long dbID) {
  Analysis *anal;
  char qStr[256];
  MYSQL_RES *results;
  MYSQL_ROW row;

  if (IDHash_contains(aa->analCache, dbID)) {
    anal = (Analysis *)IDHash_getValue(aa->analCache, dbID);

  } else {
    sprintf(qStr,
      "SELECT analysis_id, logic_name,"
      "       program, program_version, program_file,"
      "       db, db_version, db_file,"
      "       module, module_version,"
      "       gff_source, gff_feature,"
      "       created, parameters"
      " FROM   analysis"
      " WHERE  analysis_id = %d", dbID);
  
    results = aa->prepare((BaseAdaptor *)aa,qStr,strlen(qStr));
  
    row = mysql_fetch_row(results);
    if( row == NULL ) {
      return NULL;
    }
  
    anal = AnalysisAdaptor_analysisFromRow(aa, row);

    IDHash_add(aa->analCache, Analysis_getDbID(anal), anal);
    StringHash_add(aa->logicNameCache, Analysis_getLogicName(anal), anal);
  }

  return anal;
}

Analysis *AnalysisAdaptor_fetchByLogicName(AnalysisAdaptor *aa, char *logicName) {
  Analysis *anal;
  char qStr[256];
  MYSQL_RES *results;
  MYSQL_ROW row;

  if (StringHash_contains(aa->logicNameCache, logicName)) {
    anal = (Analysis *)StringHash_getValue(aa->logicNameCache, logicName);
    printf("Using logicName cache\n");

  } else {
    sprintf(qStr,
      "SELECT analysis_id, logic_name,"
      "       program, program_version, program_file,"
      "       db, db_version, db_file,"
      "       module, module_version,"
      "       gff_source, gff_feature,"
      "       created, parameters"
      " FROM   analysis"
      " WHERE  logic_name = '%s'", logicName);
  
    results = aa->prepare((BaseAdaptor *)aa,qStr,strlen(qStr));
  
    row = mysql_fetch_row(results);
    if( row == NULL ) {
      return NULL;
    }
  
    anal = AnalysisAdaptor_analysisFromRow(aa, row);

    IDHash_add(aa->analCache, Analysis_getDbID(anal), anal);
    StringHash_add(aa->logicNameCache, Analysis_getLogicName(anal), anal);
  }

  return anal;
}

Analysis *AnalysisAdaptor_analysisFromRow(AnalysisAdaptor *aa, MYSQL_ROW row) {
  Analysis *anal = Analysis_new();

  Analysis_setAdaptor(anal,(BaseAdaptor *)aa);
  Analysis_setDbID(anal,MysqlUtil_getLong(row,0));
  if (row[1]) Analysis_setLogicName(anal,MysqlUtil_getString(row,1));
  if (row[2]) Analysis_setProgram(anal,MysqlUtil_getString(row,2));
  if (row[3]) Analysis_setProgramVersion(anal,MysqlUtil_getInt(row,3));
  if (row[4]) Analysis_setProgramFile(anal,MysqlUtil_getString(row,4));
  if (row[5]) Analysis_setDb(anal,MysqlUtil_getString(row,5));
  if (row[6]) Analysis_setDbVersion(anal,MysqlUtil_getInt(row,6));
  if (row[7]) Analysis_setDbFile(anal,MysqlUtil_getString(row,7));
  if (row[8]) Analysis_setModule(anal,MysqlUtil_getString(row,8));
  if (row[9]) Analysis_setModuleVersion(anal,MysqlUtil_getInt(row,9));
  if (row[10]) Analysis_setGFFSource(anal,MysqlUtil_getString(row,10));
  if (row[11]) Analysis_setGFFFeature(anal,MysqlUtil_getString(row,11));
  if (row[12]) Analysis_setCreated(anal,MysqlUtil_getString(row,12));
  if (row[13]) Analysis_setParameters(anal,MysqlUtil_getString(row,13));

  return anal;
}
