#include "SequenceAdaptor.h"
#include "BaseAdaptor.h"
#include "MysqlUtil.h"
#include "SeqUtil.h"
#include "RawContig.h"


SequenceAdaptor *SequenceAdaptor_new(DBAdaptor *dba) {
  SequenceAdaptor *sa;

  if ((sa = (SequenceAdaptor *)calloc(1,sizeof(SequenceAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for SequenceAdaptor\n");
    return NULL;
  }
  BaseAdaptor_init((BaseAdaptor *)sa, dba, SEQUENCE_ADAPTOR);

  return sa;
}

char *SequenceAdaptor_fetchByRawContigStartEndStrand(SequenceAdaptor *sa, 
                                                   RawContig *rc,
                                                   int start,
                                                   int end,
                                                   char strand) {
  char qStr[256];
  MYSQL_RES *results;
  MYSQL_ROW row;

  if( start < 1 ) {
    fprintf(stderr,"ERROR: Wrong parameters to SequenceAdaptor fetch\n");
    return NULL;
  }

  if( end == -1 ) {

    sprintf(qStr,
            "SELECT c.length, SUBSTRING( d.sequence, %d )"
              " FROM dna d, contig c"
              " WHERE d.dna_id = c.dna_id"
              "  AND c.contig_id = %d", start, RawContig_getDbID(rc));

  } else {
    int length = end - start + 1;
    if( length < 1 ) {
      fprintf(stderr,"ERROR: Wrong parameters to SequenceAdaptor fetch\n");
      return NULL;
    }

    sprintf(qStr,"SELECT c.length,"
                 " SUBSTRING( d.sequence, %d, %d )"
                 " FROM dna d, contig c"
                 " WHERE d.dna_id = c.dna_id"
                 " AND c.contig_id = %d",start,length,RawContig_getDbID(rc));
  }

  if(  DBAdaptor_getDNADBAdaptor(sa->dba) ) {
    results = DBAdaptor_prepare(DBAdaptor_getDNADBAdaptor(sa->dba),qStr,
                                                          strlen(qStr));
  } else {
    results = sa->prepare( (BaseAdaptor *)sa, qStr, strlen(qStr) );
  }

  row = mysql_fetch_row(results);

  if( row != NULL) {
    int length   = MysqlUtil_getInt(row,0);
    char *seqStr = MysqlUtil_getString(row,1);

    /* Is this necessary ????? $seq =~ s/\s//g; */
    if( strand == -1 ) {
      return SeqUtil_reverseComplement( seqStr, strlen(seqStr) );
    } else {
      return seqStr;
    }
  } else {
    return NULL;
  }
}
