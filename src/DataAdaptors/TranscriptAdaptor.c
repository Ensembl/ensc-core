#include "TranscriptAdaptor.h"
#include "ExonAdaptor.h"
#include "BaseAdaptor.h"
#include "MysqlUtil.h"


TranscriptAdaptor *TranscriptAdaptor_new(DBAdaptor *dba) {
  TranscriptAdaptor *ta;

  if ((ta = (TranscriptAdaptor *)calloc(1,sizeof(TranscriptAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for TranscriptAdaptor\n");
    return NULL;
  }
  BaseAdaptor_init((BaseAdaptor *)ta, dba, TRANSCRIPT_ADAPTOR);

  return ta;
}

Transcript *TranscriptAdaptor_fetchByDbID(TranscriptAdaptor *ta, long dbID) {
  Transcript *trans;
  char qStr[256];
  MYSQL_RES *results;
  MYSQL_ROW row;
  ExonAdaptor *ea;
  int seen = 0;

  trans = Transcript_new();
  ea = DBAdaptor_getExonAdaptor(ta->dba);

  sprintf(qStr,
    "SELECT exon_id"
    " FROM   exon_transcript"
    " WHERE  transcript_id = %d"
    " ORDER BY rank",dbID);

  results = ta->prepare((BaseAdaptor *)ta,qStr,strlen(qStr));

  while(row = mysql_fetch_row(results)) {
    Exon *exon = ExonAdaptor_fetchByDbID(ea, MysqlUtil_getLong(row,0));
    Transcript_addExon(trans,exon);
    seen = 1;
  }

  if (!seen) {
    fprintf(stderr,"ERROR: transcript %d is not present in db",dbID);
    return NULL;
  }
  Transcript_setDbID(trans, dbID);
  Transcript_setAdaptor(trans, (BaseAdaptor *)ta );

  sprintf(qStr,
    "SELECT translation_id"
    " FROM  transcript"
    " WHERE  transcript_id = %d",dbID);

  results = ta->prepare((BaseAdaptor *)ta,qStr,strlen(qStr));

  row = mysql_fetch_row(results);
  if( row != NULL ) {
    Transcript_setTranslationId(trans, MysqlUtil_getLong(row,0));
  }

  return trans;
}

int TranscriptAdaptor_getStableEntryInfo(TranscriptAdaptor *ta, Transcript *transcript) {
  char qStr[256];
  MYSQL_RES *results;
  MYSQL_ROW row;

  if( !transcript ) {
    fprintf(stderr, "ERROR: TranscriptAdaptor_getStableEntryInfo needs a transcript object\n");
    exit(1);
  }

  sprintf(qStr,
          "SELECT stable_id, version"
          " FROM transcript_stable_id"
          " WHERE transcript_id = %d",Transcript_getDbID(transcript));

  results = ta->prepare((BaseAdaptor *)ta,qStr,strlen(qStr));

  row = mysql_fetch_row(results);
  if( row == NULL ) {
    fprintf(stderr,"WARNING: Failed fetching stable id info\n");
    return 0;
  }

  Transcript_setStableId(transcript,MysqlUtil_getString(row,0));
  Transcript_setVersion(transcript,MysqlUtil_getInt(row,1));

  return 1;
}

