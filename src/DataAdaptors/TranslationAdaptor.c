#include "TranslationAdaptor.h"
#include "BaseAdaptor.h"
#include "MysqlUtil.h"
#include "Exon.h"


TranslationAdaptor *TranslationAdaptor_new(DBAdaptor *dba) {
  TranslationAdaptor *ta;

  if ((ta = (TranslationAdaptor *)calloc(1,sizeof(TranslationAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for TranslationAdaptor\n");
    return NULL;
  }
  BaseAdaptor_init((BaseAdaptor *)ta, dba, TRANSLATION_ADAPTOR);

  return ta;
}

Translation *TranslationAdaptor_fetchByDbID(TranslationAdaptor *ta, long dbID, Transcript *transcript) {
  Translation *translation;
  char qStr[256];
  MYSQL_RES *results;
  MYSQL_ROW row;
  Exon *startExon = NULL;
  Exon *endExon = NULL;
  int i;
  long startExonId;
  long endExonId;

  if (!transcript) {
    fprintf(stderr, "ERROR: Translations make no sense outside of their " 
                    "parent Transcript objects. You must retrieve "
                    "with Transcript parent");
    exit(1);
  }

  sprintf(qStr,
    "SELECT translation_id tlid, seq_start,"
    "       start_exon_id, seq_end, end_exon_id"
    " FROM   translation"
    " WHERE  translation_id = %d", dbID);

  results = ta->prepare((BaseAdaptor *)ta,qStr,strlen(qStr));

  row = mysql_fetch_row(results);
  if( row == NULL ) {
    return NULL;
  }

  translation = Translation_new();

  Translation_setStart(translation, MysqlUtil_getInt(row,1));
  startExonId = MysqlUtil_getLong(row,2);
  Translation_setEnd(translation, MysqlUtil_getInt(row,3));
  endExonId = MysqlUtil_getLong(row,4);

  for (i=0; i<Transcript_getExonCount(transcript); i++) {
    Exon *exon = Transcript_getExonAt(transcript,i);
    if(Exon_getDbID(exon) == startExonId) {
      startExon = exon;
    }

    if(Exon_getDbID(exon) == endExonId) {
      endExon = exon;
    }
  }
  if (!startExon || !endExon) {
    fprintf(stderr,"ERROR: Could not find start or end exon in transcript\n");
    exit(1);
  }

  Translation_setStartExon(translation, startExon);
  Translation_setEndExon(translation, endExon);
  Translation_setDbID(translation, dbID);
  Translation_setAdaptor(translation, (BaseAdaptor *)ta);

  return translation;
}

int TranslationAdaptor_getStableEntryInfo(TranslationAdaptor *ta, Translation *translation) {
  char qStr[256];
  MYSQL_RES *results;
  MYSQL_ROW row;

  if( !translation ) {
    fprintf(stderr, "ERROR: TranslationAdaptor_getStableEntryInfo needs a translation object\n");
    exit(1);
  }

  sprintf(qStr,
          "SELECT stable_id, version"
          " FROM translation_stable_id"
          " WHERE translation_id = %d",Translation_getDbID(translation));

  results = ta->prepare((BaseAdaptor *)ta,qStr,strlen(qStr));

  row = mysql_fetch_row(results);
  if( row == NULL ) {
    fprintf(stderr,"WARNING: Failed fetching stable id info\n");
    return 0;
  }

  Translation_setStableId(translation,MysqlUtil_getString(row,0));
  Translation_setVersion(translation,MysqlUtil_getInt(row,1));

  return 1;
}

