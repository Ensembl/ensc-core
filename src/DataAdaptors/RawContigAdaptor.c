#include "RawContigAdaptor.h"
#include "BaseAdaptor.h"
#include "MysqlUtil.h"
#include "IDHash.h"


RawContigAdaptor *RawContigAdaptor_new(DBAdaptor *dba) {
  RawContigAdaptor *rca;

  if ((rca = (RawContigAdaptor *)calloc(1,sizeof(RawContigAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for RawContigAdaptor\n");
    return NULL;
  }
  BaseAdaptor_init((BaseAdaptor *)rca, dba, RAWCONTIG_ADAPTOR);

  rca->rawContigCache = IDHash_new(IDHASH_SMALL);

  return rca;
}

/* Lazy filling */
RawContig *RawContigAdaptor_fetchByDbID(RawContigAdaptor *rca, long dbID) {
  RawContig *rawContig;

  if (IDHash_contains(rca->rawContigCache,dbID)) {
    rawContig = (RawContig *)IDHash_getValue(rca->rawContigCache, dbID);
  } else {
    rawContig = RawContig_new();
    RawContig_setDbID(rawContig,dbID);
    IDHash_add(rca->rawContigCache,dbID,rawContig);
  }

  return rawContig;
}

void RawContigAdaptor_fetchAttributes(RawContigAdaptor *rca, RawContig *rc) {
  char qStr[256];
  MYSQL_RES *results;
  MYSQL_ROW row;

  sprintf(qStr,
          "SELECT contig_id, name, clone_id, length,"
          " embl_offset, dna_id "
          "FROM contig "
          "WHERE contig_id = %d", RawContig_getDbID(rc));

  results = rca->prepare((BaseAdaptor *)rca,qStr,strlen(qStr));

  row = mysql_fetch_row(results);
  if( row == NULL ) {
    fprintf(stderr,"ERROR: Failed fetching contig attributes\n");
    return;
  }

  RawContigAdaptor_fillRawContigWithRow(rca, rc, row);

  return;
}


void RawContigAdaptor_fillRawContigWithRow(RawContigAdaptor *rca, RawContig *rc, MYSQL_ROW row) {

  RawContig_setName(rc,MysqlUtil_getString(row,1));
  RawContig_setCloneID(rc,MysqlUtil_getLong(row,2));
  RawContig_setLength(rc,MysqlUtil_getInt(row,3));
  RawContig_setEMBLOffset(rc,MysqlUtil_getInt(row,4));

  return;
}
