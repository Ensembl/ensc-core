#include "AnalysisAdaptor.h"
#include "BaseAdaptor.h"
#include "MysqlUtil.h"
#include "StatementHandle.h"
#include "ResultRow.h"


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
  StatementHandle *sth;
  ResultRow *row;

  sprintf(qStr,
    "SELECT analysis_id, logic_name,"
    "       program, program_version, program_file,"
    "       db, db_version, db_file,"
    "       module, module_version,"
    "       gff_source, gff_feature,"
    "       created, parameters"
    " FROM   analysis");

  sth = aa->prepare((BaseAdaptor *)aa,qStr,strlen(qStr));
  sth->execute(sth);

  while (row = sth->fetchRow(sth)) {
    Analysis *anal = AnalysisAdaptor_analysisFromRow(aa, row);

    IDHash_add(aa->analCache, Analysis_getDbID(anal), anal);
    StringHash_add(aa->logicNameCache, Analysis_getLogicName(anal), anal);
  }

  sth->finish(sth);

  return (Analysis **)IDHash_getValues(aa->analCache);
}


Analysis *AnalysisAdaptor_fetchByDbID(AnalysisAdaptor *aa, int64 dbID) {
  Analysis *anal;
  char qStr[256];
  StatementHandle *sth;
  ResultRow *row;

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
      " WHERE  analysis_id = "
      INT64FMTSTR, dbID);
  
    
    sth = aa->prepare((BaseAdaptor *)aa,qStr,strlen(qStr));
  
    sth->execute(sth);

    row = sth->fetchRow(sth);
    if( row == NULL ) {
      sth->finish(sth);
      return NULL;
    }
  
    anal = AnalysisAdaptor_analysisFromRow(aa, row);
    sth->finish(sth);

    IDHash_add(aa->analCache, Analysis_getDbID(anal), anal);
    StringHash_add(aa->logicNameCache, Analysis_getLogicName(anal), anal);
  }

  return anal;
}

Analysis *AnalysisAdaptor_fetchByLogicName(AnalysisAdaptor *aa, char *logicName) {
  Analysis *anal;
  char qStr[256];
  StatementHandle *sth;
  ResultRow *row;

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
  
    sth = aa->prepare((BaseAdaptor *)aa,qStr,strlen(qStr));
    sth->execute(sth);
  
    row = sth->fetchRow(sth);
    if( row == NULL ) {
      sth->finish(sth);
      return NULL;
    }
  
    anal = AnalysisAdaptor_analysisFromRow(aa, row);
    sth->finish(sth);

    IDHash_add(aa->analCache, Analysis_getDbID(anal), anal);
    StringHash_add(aa->logicNameCache, Analysis_getLogicName(anal), anal);
  }

  return anal;
}

Analysis *AnalysisAdaptor_analysisFromRow(AnalysisAdaptor *aa, ResultRow *row) {
  Analysis *anal = Analysis_new();

  Analysis_setAdaptor(anal,(BaseAdaptor *)aa);
  Analysis_setDbID(anal,row->getLongLongAt(row,0));
  if (row->col(row,1))       Analysis_setLogicName (anal, row->getStringAt(row,1));
  if (row->col(row,2))         Analysis_setProgram (anal, row->getStringAt(row,2));
  if (row->col(row,3))  Analysis_setProgramVersion (anal, row->getIntAt(row,3));
  if (row->col(row,4))     Analysis_setProgramFile (anal, row->getStringAt(row,4));
  if (row->col(row,5))              Analysis_setDb (anal, row->getStringAt(row,5));
  if (row->col(row,6))       Analysis_setDbVersion (anal, row->getIntAt(row,6));
  if (row->col(row,7))          Analysis_setDbFile (anal, row->getStringAt(row,7));
  if (row->col(row,8))          Analysis_setModule (anal, row->getStringAt(row,8));
  if (row->col(row,9))   Analysis_setModuleVersion (anal, row->getIntAt(row,9));
  if (row->col(row,10))      Analysis_setGFFSource (anal, row->getStringAt(row,10));
  if (row->col(row,11))     Analysis_setGFFFeature (anal, row->getStringAt(row,11));
  if (row->col(row,12))        Analysis_setCreated (anal, row->getStringAt(row,12));
  if (row->col(row,13))     Analysis_setParameters (anal, row->getStringAt(row,13));

  return anal;
}
