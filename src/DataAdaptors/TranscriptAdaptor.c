#include "TranscriptAdaptor.h"
#include "ExonAdaptor.h"
#include "BaseAdaptor.h"
#include "MysqlUtil.h"

#include "StatementHandle.h"
#include "ResultRow.h"

TranscriptAdaptor *TranscriptAdaptor_new(DBAdaptor *dba) {
  TranscriptAdaptor *ta;

  if ((ta = (TranscriptAdaptor *)calloc(1,sizeof(TranscriptAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for TranscriptAdaptor\n");
    return NULL;
  }
  BaseAdaptor_init((BaseAdaptor *)ta, dba, TRANSCRIPT_ADAPTOR);

  return ta;
}

Transcript *TranscriptAdaptor_fetchByDbID(TranscriptAdaptor *ta, IDType dbID) {
  Transcript *trans;
  char qStr[256];
  StatementHandle *sth;
  ResultRow *row;
  ExonAdaptor *ea;
  int seen = 0;

  trans = Transcript_new();
  ea = DBAdaptor_getExonAdaptor(ta->dba);

  sprintf(qStr,
    "SELECT exon_id"
    " FROM   exon_transcript"
    " WHERE  transcript_id = "
    IDFMTSTR
    " ORDER BY rank",dbID);

  sth = ta->prepare((BaseAdaptor *)ta,qStr,strlen(qStr));
  sth->execute(sth);

  while(row = sth->fetchRow(sth)) {
    Exon *exon = ExonAdaptor_fetchByDbID(ea, row->getLongLongAt(row,0));
    Transcript_addExon(trans,exon);
    seen = 1;
  }
  sth->finish(sth);

  if (!seen) {
    fprintf(stderr,"ERROR: transcript " IDFMTSTR "is not present in db",dbID);
    return NULL;
  }
  Transcript_setDbID(trans, dbID);
  Transcript_setAdaptor(trans, (BaseAdaptor *)ta );

  sprintf(qStr,
    "SELECT translation_id"
    " FROM  transcript"
    " WHERE  transcript_id = " IDFMTSTR,dbID);

  sth = ta->prepare((BaseAdaptor *)ta,qStr,strlen(qStr));
  sth->execute(sth);

  row = sth->fetchRow(sth);
  if( row != NULL ) {
    Transcript_setTranslationId(trans, row->getLongLongAt(row,0));
  }
  sth->finish(sth);

  return trans;
}

int TranscriptAdaptor_getStableEntryInfo(TranscriptAdaptor *ta, Transcript *transcript) {
  char qStr[256];
  StatementHandle *sth;
  ResultRow *row;

  if( !transcript ) {
    fprintf(stderr, "ERROR: TranscriptAdaptor_getStableEntryInfo needs a transcript object\n");
    exit(1);
  }

  sprintf(qStr,
          "SELECT stable_id, version"
          " FROM transcript_stable_id"
          " WHERE transcript_id = " IDFMTSTR, Transcript_getDbID(transcript));

  sth = ta->prepare((BaseAdaptor *)ta,qStr,strlen(qStr));
  sth->execute(sth);

  row = sth->fetchRow(sth);
  if( row == NULL ) {
    fprintf(stderr,"WARNING: Failed fetching stable id info\n");
    sth->finish(sth);
    return 0;
  }

  Transcript_setStableId(transcript,row->getStringAt(row,0));
  Transcript_setVersion(transcript,row->getIntAt(row,1));

  sth->finish(sth);

  return 1;
}

