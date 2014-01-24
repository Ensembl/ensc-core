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
  
// Just populating caches, so free returned array
  Analysis **all = AnalysisAdaptor_fetchAll(aa);
  free(all);

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

  while ((row = sth->fetchRow(sth))) {
    Analysis *anal = AnalysisAdaptor_analysisFromRow(aa, row);

    IDHash_add(aa->analCache, Analysis_getDbID(anal), anal);
    StringHash_add(aa->logicNameCache, Analysis_getLogicName(anal), anal);
  }

  sth->finish(sth);

  return (Analysis **)IDHash_getValues(aa->analCache);
}


Analysis *AnalysisAdaptor_fetchByDbID(AnalysisAdaptor *aa, IDType dbID) {
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
      IDFMTSTR, dbID);
  
    
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
  char qStr[1024];
  StatementHandle *sth;
  ResultRow *row;

  if (StringHash_contains(aa->logicNameCache, logicName)) {
    anal = (Analysis *)StringHash_getValue(aa->logicNameCache, logicName);

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
  if (row->col(row,1))       Analysis_setLogicName (anal, row->getStringCopyAt(row,1));
  if (row->col(row,2))         Analysis_setProgram (anal, row->getStringCopyAt(row,2));
  if (row->col(row,3))  Analysis_setProgramVersion (anal, row->getIntAt(row,3));
  if (row->col(row,4))     Analysis_setProgramFile (anal, row->getStringCopyAt(row,4));
  if (row->col(row,5))              Analysis_setDb (anal, row->getStringCopyAt(row,5));
  if (row->col(row,6))       Analysis_setDbVersion (anal, row->getIntAt(row,6));
  if (row->col(row,7))          Analysis_setDbFile (anal, row->getStringCopyAt(row,7));
  if (row->col(row,8))          Analysis_setModule (anal, row->getStringCopyAt(row,8));
  if (row->col(row,9))   Analysis_setModuleVersion (anal, row->getIntAt(row,9));
  if (row->col(row,10))      Analysis_setGFFSource (anal, row->getStringCopyAt(row,10));
  if (row->col(row,11))     Analysis_setGFFFeature (anal, row->getStringCopyAt(row,11));
  if (row->col(row,12))        Analysis_setCreated (anal, row->getStringCopyAt(row,12));
  if (row->col(row,13))     Analysis_setParameters (anal, row->getStringCopyAt(row,13));

  return anal;
}

IDType AnalysisAdaptor_store(AnalysisAdaptor *aa, Analysis *analysis) {
  IDType dbID = Analysis_getDbID(analysis);
  StatementHandle *sth;
  char qStr[4096];
  char fmtStr[1024];
  
  //fprintf(stderr,"Analysis dbId = " IDFMTSTR "\n", dbID);

  if (dbID && (AnalysisAdaptor *)Analysis_getAdaptor(analysis) == aa) {
    //fprintf(stderr, "Returning - think have adaptor and dbID\n");
    return dbID;
  }


  dbID = 0;

  if ((dbID = AnalysisAdaptor_analysisExists(aa, analysis))) {
    //fprintf(stderr, "Returning - think exists\n");
    Analysis_setAdaptor(analysis,(BaseAdaptor *)aa);
    Analysis_setDbID(analysis,dbID);
    return dbID;
  }
 
  if (Analysis_getLogicName(analysis) == NULL) {
    fprintf(stderr,"Must have a logic name on the analysis object");
    exit(1);
  }


  sprintf(fmtStr,
    "INSERT IGNORE INTO analysis"
    " SET created = %%s,"
    "     logic_name = %%s,"
    "     db = %%s,"
    "     db_version = %%s,"
    "     db_file = %%s,"
    "     program = %%s,"
    "     program_version = %%s,"
    "     program_file = %%s,"
    "     parameters = %%s,"
    "     module = %%s,"
    "     module_version = %%s,"
    "     gff_source = %%s,"
    "     gff_feature = %%s");

  char createdQStr[1024];
  Analysis_getCreated(analysis) ? sprintf(createdQStr,"'%s'", Analysis_getCreated(analysis)) : sprintf(createdQStr, "now()");

  char logicNameQStr[1024];
  Analysis_getLogicName(analysis) ? sprintf(logicNameQStr,"'%s'", Analysis_getLogicName(analysis)) : sprintf(logicNameQStr, "NULL");

  char dbQStr[1024];
  Analysis_getDb(analysis) ? sprintf(dbQStr,"'%s'", Analysis_getDb(analysis)) : sprintf(dbQStr, "NULL");

  char dbVerQStr[1024];
  Analysis_getDbVersion(analysis) ? sprintf(dbVerQStr,"%d", Analysis_getDbVersion(analysis)) : sprintf(dbVerQStr, "NULL");

  char dbFileQStr[1024];
  Analysis_getDbFile(analysis) ? sprintf(dbFileQStr,"'%s'", Analysis_getDbFile(analysis)) : sprintf(dbFileQStr, "NULL");

  char programQStr[1024];
  Analysis_getProgram(analysis) ? sprintf(programQStr,"'%s'", Analysis_getProgram(analysis)) : sprintf(programQStr, "NULL");

  char progVerQStr[1024];
  Analysis_getProgramVersion(analysis) ? sprintf(progVerQStr,"%d", Analysis_getProgramVersion(analysis)) : sprintf(progVerQStr, "NULL");

  char progFileQStr[1024];
  Analysis_getProgramFile(analysis) ? sprintf(progFileQStr,"'%s'", Analysis_getProgramFile(analysis)) : sprintf(progFileQStr, "NULL");

  char paramQStr[1024];
  Analysis_getParameters(analysis) ? sprintf(paramQStr,"'%s'", Analysis_getParameters(analysis)) : sprintf(paramQStr, "NULL");

  char moduleQStr[1024];
  Analysis_getModule(analysis) ? sprintf(moduleQStr,"'%s'", Analysis_getModule(analysis)) : sprintf(moduleQStr, "NULL");

  char modVerQStr[1024];
  Analysis_getModuleVersion(analysis) ? sprintf(modVerQStr,"%d", Analysis_getModuleVersion(analysis)) : sprintf(modVerQStr, "NULL");

  char gffSrcQStr[1024];
  Analysis_getGFFSource(analysis) ? sprintf(gffSrcQStr,"'%s'", Analysis_getGFFSource(analysis)) : sprintf(gffSrcQStr, "NULL");
 
  char gffFeatQStr[1024];
  Analysis_getGFFFeature(analysis) ? sprintf(gffFeatQStr,"'%s'", Analysis_getGFFFeature(analysis)) : sprintf(gffFeatQStr, "NULL");

  sprintf(qStr, fmtStr, createdQStr, logicNameQStr, dbQStr, dbVerQStr, dbFileQStr, programQStr, progVerQStr, progFileQStr, paramQStr, moduleQStr, modVerQStr, gffSrcQStr, gffFeatQStr);

  sth = aa->prepare((BaseAdaptor *)aa,qStr,strlen(qStr));
  int nRowInserted = sth->execute(sth);
  dbID = sth->getInsertId(sth);
  sth->finish(sth);

  if (!nRowInserted) {
    // if we didn't insert it should hopefully mean its already there so fetch it
    Analysis *tmpAnal = AnalysisAdaptor_fetchByLogicName(aa, Analysis_getLogicName(analysis));
    if (tmpAnal == NULL) {
      fprintf(stderr,"Failed storing analysis\n");
      exit(1);
    }
    dbID = Analysis_getDbID(tmpAnal);
    Analysis_setCreated(analysis, Analysis_getCreated(tmpAnal));

  } else if (!Analysis_getCreated(analysis)) {

    if (dbID) {
      ResultRow *row;

      sprintf(qStr,"SELECT created FROM analysis WHERE analysis_id = " IDFMTSTR, dbID); 
      sth = aa->prepare((BaseAdaptor *)aa,qStr,strlen(qStr));
      sth->execute(sth);
      row = sth->fetchRow(sth);
      Analysis_setCreated(analysis, row->getStringAt(row,0));
      sth->finish(sth);
    }
  }

  //fprintf(stderr,"dbID = "IDFMTSTR"\n", dbID);

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

IDType AnalysisAdaptor_analysisExists(AnalysisAdaptor *aa, Analysis *anal) {
  IDType *keys;
  int i;
  IDType cacheId = 0;


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
    Analysis *compAnal = (Analysis *)IDHash_getValue(aa->analCache, keys[i]);

    if (Analysis_compare(compAnal, anal) == 0) {
      //fprintf(stderr,"Think analysis %p %s matches analysis %p %s\n", anal, Analysis_getLogicName(anal), compAnal, Analysis_getLogicName(compAnal));
      cacheId = keys[i];
      break;
    }
  }
  free(keys);
  return cacheId;
}
