#include "ExonAdaptor.h"
#include "BaseAdaptor.h"
#include "MysqlUtil.h"
#include "RawContigAdaptor.h"
#include "IDHash.h"

#include "StatementHandle.h"
#include "ResultRow.h"


Exon *ExonAdaptor_exonFromRow(ExonAdaptor *ea, ResultRow *row);

ExonAdaptor *ExonAdaptor_new(DBAdaptor *dba) {
  ExonAdaptor *ea;

  if ((ea = (ExonAdaptor *)calloc(1,sizeof(ExonAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for ExonAdaptor\n");
    return NULL;
  }
  BaseAdaptor_init((BaseAdaptor *)ea, dba, EXON_ADAPTOR);

  return ea;
}

Exon *ExonAdaptor_fetchByDbID(ExonAdaptor *ea, int64 dbID) {
  Exon *exon;
  char qStr[256];
  StatementHandle *sth;
  ResultRow *row;

  sprintf(qStr,
    "SELECT exon_id"
    " , contig_id"
    " , contig_start"
    " , contig_end"
    " , contig_strand"
    " , phase"
    " , end_phase"
    " , sticky_rank"
    " FROM   exon"
    " WHERE  exon_id = %d"
    " ORDER BY sticky_rank DESC", 
    dbID);

  sth = ea->prepare((BaseAdaptor *)ea,qStr,strlen(qStr));

  row = sth->fetchRow(sth);
  if( row == NULL ) {
    sth->finish(sth);
    return NULL;
  }

  exon = ExonAdaptor_exonFromResults(ea, sth, row);
  sth->finish(sth);

  return exon;
}

Exon *ExonAdaptor_exonFromResults(ExonAdaptor *ea, StatementHandle *sth, ResultRow *row) {
  Exon *exon;
  int maxRank = row->getIntAt(row,7);

  if (maxRank > 1) {
    int stickyLength = 0;
    Exon *component;
    
    fprintf(stderr, "ERROR: Sticky exons not implemented yet\n");

    // sticky exon
    exon = Exon_new();
    Exon_setDbID(exon, row->getLongLongAt(row,0));

    // make first component exon
    component = ExonAdaptor_exonFromRow(ea, row);

    Exon_addComponentExon(exon,component);
    stickyLength += Exon_getLength(component);

    Exon_setPhase(exon,Exon_getPhase(component));
    Exon_setEndPhase(exon,Exon_getEndPhase(component));
    Exon_setAdaptor(exon,(BaseAdaptor *)ea);

    // continue while loop until we hit sticky_rank 1
    while( row = sth->fetchRow(sth)) {
      component = ExonAdaptor_exonFromRow(ea, row);
  
      Exon_addComponentExon(exon,component);
      stickyLength += Exon_getLength(component);

      if( Exon_getStickyRank(component) == 1 ) {
        Exon_setContig(exon, Exon_getContig(component));
        break;
      }
    }

    Exon_sortByStickyRank(exon);

    Exon_setStart(exon,1);
    Exon_setEnd(exon,stickyLength);
    Exon_setStrand(exon, 1 );

  } else {
    exon = ExonAdaptor_exonFromRow(ea, row);
  }

  return exon;
}

Exon *ExonAdaptor_exonFromRow(ExonAdaptor *ea, ResultRow *row) {
  Exon *exon = Exon_new();
  RawContigAdaptor *rca;
  RawContig *rc;

  Exon_setDbID(exon,row->getLongLongAt(row,0));
  Exon_setStart(exon,row->getLongAt(row,2));
  Exon_setEnd(exon,row->getLongAt(row,3));
  Exon_setStrand(exon,row->getIntAt(row,4));
  Exon_setPhase(exon,row->getIntAt(row,5));
  Exon_setEndPhase(exon,row->getIntAt(row,6));
  Exon_setStickyRank(exon,row->getIntAt(row,7));

  Exon_setAdaptor(exon,(BaseAdaptor *)ea);

  rca = DBAdaptor_getRawContigAdaptor(ea->dba);
  rc = RawContigAdaptor_fetchByDbID(rca,row->getLongLongAt(row,1));

  Exon_setContig(exon,rc);

  return exon; 
}

int ExonAdaptor_fetchAllByGeneId(ExonAdaptor *ea, int64 geneId, Exon ***retExons) {
  Exon **exons;
  char qStr[512];
  StatementHandle *sth;
  ResultRow *row;
  IDHash *exonHash = IDHash_new(IDHASH_SMALL);
  int nExon;

  if( !geneId ) {
    fprintf(stderr,"ERROR: Gene dbID not defined\n");
  }

  sprintf(qStr,
    "SELECT  e.exon_id"
    "  , e.contig_id"
    "  , e.contig_start"
    "  , e.contig_end"
    "  , e.contig_strand"
    "  , e.phase"
    "  , e.end_phase"
    "  , e.sticky_rank"
    " FROM exon e"
    "  , exon_transcript et"
    "  , transcript t"
    " WHERE t.gene_id = %d"
    "  AND et.transcript_id = t.transcript_id"
    "  AND e.exon_id = et.exon_id"
    " ORDER BY t.transcript_id,e.exon_id"
    "  , e.sticky_rank DESC",geneId);

  sth = ea->prepare((BaseAdaptor *)ea, qStr, strlen(qStr));
  sth->execute(sth);

  while( row = sth->fetchRow(sth)) {
    if( ! IDHash_contains(exonHash,row->getLongLongAt(row,0))) {
      Exon *exon = ExonAdaptor_exonFromResults(ea,sth,row);

      IDHash_add(exonHash,Exon_getDbID(exon),exon);
    }
  }
  sth->finish(sth);

  *retExons = (Exon **)IDHash_getValues(exonHash);

  nExon = IDHash_getNumValues(exonHash);

  IDHash_free(exonHash,NULL);

  return nExon;

}

int ExonAdaptor_getStableEntryInfo(ExonAdaptor *ea, Exon *exon) {
  char qStr[256];
  StatementHandle *sth;
  ResultRow *row;

  if( !exon ) {
    fprintf(stderr, "ERROR: ExonAdaptor_getStableEntryInfo needs a exon object\n");
    exit(1);
  }

  sprintf(qStr,
          "SELECT stable_id, UNIX_TIMESTAMP(created),"
          "                  UNIX_TIMESTAMP(modified), version"
          " FROM exon_stable_id"
          " WHERE exon_id = %d",Exon_getDbID(exon));

  sth = ea->prepare((BaseAdaptor *)ea,qStr,strlen(qStr));

  row = sth->fetchRow(sth);
  if( row == NULL ) {
    sth->finish(sth);
    fprintf(stderr,"WARNING: Failed fetching stable id info\n");
    return 0;
  }

  Exon_setStableId(exon,row->getStringAt(row,0));
  Exon_setCreated(exon,row->getIntAt(row,1));
  Exon_setModified(exon,row->getIntAt(row,2));
  Exon_setVersion(exon,row->getIntAt(row,3));

  sth->finish(sth);

  return 1;
}

