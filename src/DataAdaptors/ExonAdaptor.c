#include "ExonAdaptor.h"
#include "BaseAdaptor.h"
#include "MysqlUtil.h"
#include "RawContigAdaptor.h"
#include "IDHash.h"


Exon *ExonAdaptor_exonFromRow(ExonAdaptor *ea, MYSQL_ROW row);

ExonAdaptor *ExonAdaptor_new(DBAdaptor *dba) {
  ExonAdaptor *ea;

  if ((ea = (ExonAdaptor *)calloc(1,sizeof(ExonAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for ExonAdaptor\n");
    return NULL;
  }
  BaseAdaptor_init((BaseAdaptor *)ea, dba, EXON_ADAPTOR);

  return ea;
}

Exon *ExonAdaptor_fetchByDbID(ExonAdaptor *ea, long dbID) {
  Exon *exon;
  char qStr[256];
  MYSQL_RES *results;
  MYSQL_ROW row;

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

  results = ea->prepare((BaseAdaptor *)ea,qStr,strlen(qStr));

  row = mysql_fetch_row(results);
  if( row == NULL ) {
    return NULL;
  }

  exon = ExonAdaptor_exonFromResults(ea, results, row);

  return exon;
}

Exon *ExonAdaptor_exonFromResults(ExonAdaptor *ea, MYSQL_RES *results, MYSQL_ROW row) {
  Exon *exon;
  int maxRank = MysqlUtil_getInt(row,7);

  if (maxRank > 1) {
    int stickyLength = 0;
    Exon *component;
    
    fprintf(stderr, "ERROR: Sticky exons not implemented yet\n");

    // sticky exon
    exon = Exon_new();
    Exon_setDbID(exon, MysqlUtil_getLong(row,0));

    // make first component exon
    component = ExonAdaptor_exonFromRow(ea, row);

    Exon_addComponentExon(exon,component);
    stickyLength += Exon_getLength(component);

    Exon_setPhase(exon,Exon_getPhase(component));
    Exon_setEndPhase(exon,Exon_getEndPhase(component));
    Exon_setAdaptor(exon,(BaseAdaptor *)ea);

    // continue while loop until we hit sticky_rank 1
    while( row = mysql_fetch_row(results)) {
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

Exon *ExonAdaptor_exonFromRow(ExonAdaptor *ea, MYSQL_ROW row) {
  Exon *exon = Exon_new();
  RawContigAdaptor *rca;
  RawContig *rc;

  Exon_setDbID(exon,MysqlUtil_getLong(row,0));
  Exon_setStart(exon,MysqlUtil_getLong(row,2));
  Exon_setEnd(exon,MysqlUtil_getLong(row,3));
  Exon_setStrand(exon,MysqlUtil_getInt(row,4));
  Exon_setPhase(exon,MysqlUtil_getInt(row,5));
  Exon_setEndPhase(exon,MysqlUtil_getInt(row,6));
  Exon_setStickyRank(exon,MysqlUtil_getInt(row,7));
  Exon_setAdaptor(exon,(BaseAdaptor *)ea);

  rca = DBAdaptor_getRawContigAdaptor(ea->dba);
  rc = RawContigAdaptor_fetchByDbID(rca,MysqlUtil_getLong(row,1));

  Exon_setContig(exon,rc);

  return exon; 
}

int ExonAdaptor_fetchAllByGeneId(ExonAdaptor *ea, long geneId, Exon ***retExons) {
  Exon **exons;
  char qStr[512];
  MYSQL_RES *results;
  MYSQL_ROW row;
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

  results = ea->prepare((BaseAdaptor *)ea, qStr, strlen(qStr));

  while( row = mysql_fetch_row(results)) {
    if( ! IDHash_contains(exonHash,MysqlUtil_getLong(row,0))) {
      Exon *exon = ExonAdaptor_exonFromResults(ea,results,row);

      IDHash_add(exonHash,Exon_getDbID(exon),exon);
    }
  }

  *retExons = (Exon **)IDHash_getValues(exonHash);

  nExon = IDHash_getNumValues(exonHash);

  IDHash_free(exonHash,NULL);

  return nExon;

}

int ExonAdaptor_getStableEntryInfo(ExonAdaptor *ea, Exon *exon) {
  char qStr[256];
  MYSQL_RES *results;
  MYSQL_ROW row;

  if( !exon ) {
    fprintf(stderr, "ERROR: ExonAdaptor_getStableEntryInfo needs a exon object\n");
    exit(1);
  }

  sprintf(qStr,
          "SELECT stable_id, UNIX_TIMESTAMP(created),"
          "                  UNIX_TIMESTAMP(modified), version"
          " FROM exon_stable_id"
          " WHERE exon_id = %d",Exon_getDbID(exon));

  results = ea->prepare((BaseAdaptor *)ea,qStr,strlen(qStr));

  row = mysql_fetch_row(results);
  if( row == NULL ) {
    fprintf(stderr,"WARNING: Failed fetching stable id info\n");
    return 0;
  }

  Exon_setStableId(exon,MysqlUtil_getString(row,0));
  Exon_setCreated(exon,MysqlUtil_getInt(row,1));
  Exon_setModified(exon,MysqlUtil_getInt(row,2));
  Exon_setVersion(exon,MysqlUtil_getInt(row,3));

  return 1;
}

