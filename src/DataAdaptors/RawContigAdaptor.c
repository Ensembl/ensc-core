#include "RawContigAdaptor.h"
#include "BaseAdaptor.h"
#include "MysqlUtil.h"
#include "IDHash.h"

#include "StatementHandle.h"
#include "ResultRow.h"


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
RawContig *RawContigAdaptor_fetchByDbID(RawContigAdaptor *rca, IDType dbID) {
  RawContig *rawContig;

  if (IDHash_contains(rca->rawContigCache,dbID)) {
    rawContig = (RawContig *)IDHash_getValue(rca->rawContigCache, dbID);
  } else {
    rawContig = RawContig_new();
    RawContig_setDbID(rawContig,dbID);
    RawContig_setAdaptor(rawContig,(BaseAdaptor *)rca);
    IDHash_add(rca->rawContigCache,dbID,rawContig);
  }

  return rawContig;
}

void RawContigAdaptor_fetchAttributes(RawContigAdaptor *rca, RawContig *rc) {
  char qStr[256];
  StatementHandle *sth;
  ResultRow *row;

  sprintf(qStr,
          "SELECT contig_id, name, clone_id, length,"
          " embl_offset, dna_id "
          "FROM contig "
          "WHERE contig_id = " IDFMTSTR, RawContig_getDbID(rc));

  sth = rca->prepare((BaseAdaptor *)rca,qStr,strlen(qStr));
  sth->execute(sth);

  row = sth->fetchRow(sth);
  if( row == NULL ) {
    fprintf(stderr,"ERROR: Failed fetching contig attributes\n");
    return;
  }

  RawContigAdaptor_fillRawContigWithRow(rca, rc, row);
  sth->finish(sth);

  return;
}


void RawContigAdaptor_fillRawContigWithRow(RawContigAdaptor *rca, RawContig *rc, ResultRow *row) {

//NIY Freeing of name
  RawContig_setName(rc,row->getStringAt(row,1));
  RawContig_setCloneID(rc,row->getLongLongAt(row,2));
  RawContig_setLength(rc,row->getIntAt(row,3));
  RawContig_setEMBLOffset(rc,row->getIntAt(row,4));

  return;
}

RawContig *RawContigAdaptor_fetchByName(RawContigAdaptor *rca, char *name) {
  char qStr[256];
  StatementHandle *sth;
  ResultRow *row;
  RawContig *rawContig;

  sprintf(qStr,
           "SELECT contig_id, name, clone_id, length,"
                           " embl_offset, dna_id"
                           " FROM contig"
                           " WHERE name = '%s'",name );
  sth = rca->prepare((BaseAdaptor *)rca, qStr,strlen(qStr));
  sth->execute(sth);

  if ((row = sth->fetchRow(sth))) {
    rawContig = RawContig_new();
    RawContig_setDbID(rawContig,row->getLongLongAt(row,0));
    RawContig_setAdaptor(rawContig,(BaseAdaptor *)rca);
    RawContigAdaptor_fillRawContigWithRow(rca, rawContig, row);
  } else {
    fprintf(stderr,"Couldn't find contig with name %s\n",name);
    exit(1);
  }
  sth->finish(sth);

  return rawContig;
}

