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

int64 AnalysisAdaptor_store(AnalysisAdaptor *aa, Analysis *analysis) {
  int64 dbID = Analysis_getDbID(analysis);
  StatementHandle *sth;
  char qStr[1024];
  
  if (dbID && (AnalysisAdaptor *)Analysis_getAdaptor(analysis) == aa) 
    return dbID;


  dbID = 0;

  if (dbID = AnalysisAdaptor_analysisExists(aa, analysis)) {
    Analysis_setAdaptor(analysis,(BaseAdaptor *)aa);
    Analysis_setDbID(analysis,dbID);
    return dbID;
  }
 
  if (Analysis_getLogicName(analysis)) {
    fprintf(stderr,"Must have a logic name on the analysis object");
    exit(1);
  }

  if (Analysis_getCreated(analysis)) {
    sprintf(qStr,
      "INSERT INTO analysis"
      " SET created = %%d,"
      "     logic_name = %%s,"
      "     db = %%s,"
      "     db_version = %%d,"
      "     db_file = %%s,"
      "     program = %%s,"
      "     program_version = %%d,"
      "     program_file = %%s,"
      "     parameters = %%s,"
      "     module = %%s,"
      "     module_version = %%d,"
      "     gff_source = %%s,"
      "     gff_feature = %%s");
    sth = aa->prepare((BaseAdaptor *)aa,qStr,strlen(qStr));
    sth->execute(sth,
        Analysis_getCreated(analysis),
	Analysis_getLogicName(analysis),
	Analysis_getDb(analysis),
	Analysis_getDbVersion(analysis),
	Analysis_getDbFile(analysis),
	Analysis_getProgram(analysis),
	Analysis_getProgramVersion(analysis),
	Analysis_getProgramFile(analysis),
	Analysis_getParameters(analysis),
	Analysis_getModule(analysis),
	Analysis_getModuleVersion(analysis),
	Analysis_getGFFSource(analysis),
	Analysis_getGFFFeature(analysis)
      );
    dbID = sth->getInsertId(sth);
    sth->finish(sth);
  } else {
    sprintf(qStr,
      "INSERT INTO analysis"
      " SET created = now(),"
      "     logic_name = %%s,"
      "     db = %%s,"
      "     db_version = %%d,"
      "     db_file = %%s,"
      "     program = %%s,"
      "     program_version = %%d,"
      "     program_file = %%s,"
      "     parameters = %%s,"
      "     module = %%s,"
      "     module_version = %%d,"
      "     gff_source = %%s,"
      "     gff_feature = %%s");
    sth = aa->prepare((BaseAdaptor *)aa,qStr,strlen(qStr));
    sth->execute(sth,
        Analysis_getCreated(analysis),
	Analysis_getLogicName(analysis),
	Analysis_getDb(analysis),
	Analysis_getDbVersion(analysis),
	Analysis_getDbFile(analysis),
	Analysis_getProgram(analysis),
	Analysis_getProgramVersion(analysis),
	Analysis_getProgramFile(analysis),
	Analysis_getParameters(analysis),
	Analysis_getModule(analysis),
	Analysis_getModuleVersion(analysis),
	Analysis_getGFFSource(analysis),
	Analysis_getGFFFeature(analysis)
      );
    dbID = sth->getInsertId(sth);
    sth->finish(sth);

    if (dbID) {
      ResultRow *row;

      sprintf(qStr,"SELECT created FROM analysis WHERE analysis_id = " INT64FMTSTR, dbID); 
      sth = aa->prepare((BaseAdaptor *)aa,qStr,strlen(qStr));
      sth->execute(sth);
      row = sth->fetchRow(sth);
      Analysis_setCreated(analysis, row->getStringAt(row,0));
      sth->finish(sth);
    }
  }
  IDHash_add(aa->analCache, dbID, analysis);
  StringHash_add(aa->logicNameCache, Analysis_getLogicName(analysis), analysis);

  Analysis_setAdaptor(analysis,(BaseAdaptor *)aa);
  Analysis_setDbID(analysis,dbID);

  return dbID;
}

/*
sub remove {
  my ($self, $analysis) = @_;
  my $dbID;
  
  if (!defined $analysis || !ref $analysis) {
    $self->throw("called remove on AnalysisAdaptor with a [$analysis]");
  }

  unless ($dbID = $self->exists($analysis)) {
    return undef;
  }

  my $res = $self->db->db_handle->do(qq{
    DELETE from analysis
    WHERE  analysis_id = $dbID
  });
}
*/

int64 AnalysisAdaptor_analysisExists(AnalysisAdaptor *aa, Analysis *anal) {
  int64 *keys;
  int i;
  int64 cacheId = 0;


  // objects with already have this adaptor are store here.
  if ((AnalysisAdaptor *)Analysis_getAdaptor(anal) == aa) {
    if (Analysis_getDbID(anal)) {
      return Analysis_getDbID(anal);
    } else {
      fprintf(stderr,"analysis does not have an analysisId");
      exit(1);
    }
  }
  
  keys = IDHash_getKeys(aa->analCache);
  for (i=0;i<IDHash_getNumValues(aa->analCache);i++) {
    if (Analysis_compare((Analysis *)IDHash_getValue(aa->analCache, keys[i]) , anal) >= 0) {
      cacheId = keys[i];
      break;
    }
  }
  free(keys);
  return cacheId;
}
