#include "ExonAdaptor.h"
#include "DBAdaptor.h"
#include "BaseAdaptor.h"
#include "MysqlUtil.h"
#include "RawContigAdaptor.h"
#include "IDHash.h"

#include "StatementHandle.h"
#include "ResultRow.h"

#include "Class.h"
#include "BaseAlignFeature.h"


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

Exon *ExonAdaptor_fetchByDbID(ExonAdaptor *ea, IDType dbID) {
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
    " WHERE  exon_id = "
    IDFMTSTR
    " ORDER BY sticky_rank DESC", 
    dbID);

  sth = ea->prepare((BaseAdaptor *)ea,qStr,strlen(qStr));
  sth->execute(sth);

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
    StickyExon *stickyExon;
    
    fprintf(stderr, "ERROR: Sticky exons not implemented yet\n");

    // sticky exon
    stickyExon = StickyExon_new();
    exon = (Exon *)stickyExon;

    StickyExon_setDbID(exon, row->getLongLongAt(row,0));

    // make first component exon
    component = ExonAdaptor_exonFromRow(ea, row);

    StickyExon_addComponentExon(stickyExon,component);
    stickyLength += Exon_getLength(component);

    StickyExon_setPhase(stickyExon,Exon_getPhase(component));
    StickyExon_setEndPhase(stickyExon,Exon_getEndPhase(component));
    StickyExon_setAdaptor(stickyExon,(BaseAdaptor *)ea);

    // continue while loop until we hit sticky_rank 1
    while ((row = sth->fetchRow(sth))) {
      component = ExonAdaptor_exonFromRow(ea, row);
  
      StickyExon_addComponentExon(stickyExon,component);
      stickyLength += Exon_getLength(component);

      if( Exon_getStickyRank(component) == 1 ) {
        StickyExon_setContig(stickyExon, Exon_getContig(component));
        break;
      }
    }

    StickyExon_sortByStickyRank(stickyExon);

    StickyExon_setStart(stickyExon,1);
    StickyExon_setEnd(stickyExon,stickyLength);
    StickyExon_setStrand(stickyExon, 1 );

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
  // printf("Set exon id to " IDFMTSTR " from " IDFMTSTR "\n",Exon_getDbID(exon),row->getLongLongAt(row,0));
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

int ExonAdaptor_fetchAllByGeneId(ExonAdaptor *ea, IDType geneId, Exon ***retExons) {
  char qStr[512];
  StatementHandle *sth;
  ResultRow *row;
  IDHash *exonHash = IDHash_new(IDHASH_SMALL);
  int nExon;

  if( !geneId ) {
    fprintf(stderr,"ERROR: Gene dbID not defined\n");
  }

  sprintf(qStr,
    "SELECT STRAIGHT_JOIN e.exon_id"
    "  , e.contig_id"
    "  , e.contig_start"
    "  , e.contig_end"
    "  , e.contig_strand"
    "  , e.phase"
    "  , e.end_phase"
    "  , e.sticky_rank"
    " FROM transcript t"
    "  , exon_transcript et"
    "  , exon e"
    " WHERE t.gene_id = "
    IDFMTSTR
    "  AND et.transcript_id = t.transcript_id"
    "  AND e.exon_id = et.exon_id"
    " ORDER BY t.transcript_id,e.exon_id"
    "  , e.sticky_rank DESC",geneId);

  sth = ea->prepare((BaseAdaptor *)ea, qStr, strlen(qStr));
  sth->execute(sth);

  while ((row = sth->fetchRow(sth))) {
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
          " WHERE exon_id = " IDFMTSTR, Exon_getDbID(exon));

  sth = ea->prepare((BaseAdaptor *)ea,qStr,strlen(qStr));
  sth->execute(sth);

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

char *protType = "protein_align_feature";
char *dnaType = "dna_align_feature";

IDType  ExonAdaptor_store(ExonAdaptor *ea, Exon *exon) {
  StatementHandle *sth;
  StatementHandle *sth2;
  char qStr[1024];
  int i;
  IDType exonId;
  DNAAlignFeatureAdaptor *dafa;
  ProteinAlignFeatureAdaptor *pafa;
  char *type;
  int nExon;
  StickyExon *stickyExon = NULL;
  

  Class_assertType(CLASS_EXON, exon->objectType);

  if (Exon_getDbID(exon) && Exon_getAdaptor(exon) && Exon_getAdaptor(exon) == (BaseAdaptor *)ea) {
    return Exon_getDbID(exon);
  }

  if( ! Exon_getStart(exon)  || ! Exon_getEnd(exon) ||
      ! Exon_getStrand(exon) || ! Exon_getPhase(exon)) {
    fprintf(stderr,"ERROR: Exon does not have all attributes to store");
    exit(1);
  }

  // trap contig_id separately as it is likely to be a common mistake

// HACK modified duplicate of query for first exon (with no setting of exonId
  sprintf(qStr,
    "INSERT into exon (exon_id, contig_id, contig_start,"
                      "contig_end, contig_strand, phase,"
                      "end_phase, sticky_rank) "
    " VALUES ( %" IDFMTSTR ", %" IDFMTSTR ", %%d, %%d, %%d, %%d, %%d, %%d )");


  sth = ea->prepare((BaseAdaptor *)ea,qStr,strlen(qStr));

  sprintf(qStr,
    "INSERT into exon (contig_id, contig_start,"
                      "contig_end, contig_strand, phase,"
                      "end_phase, sticky_rank) "
    " VALUES ( %" IDFMTSTR ", %%d, %%d, %%d, %%d, %%d, %%d )");

  sth2 = ea->prepare((BaseAdaptor *)ea,qStr,strlen(qStr));

  exonId = 0;

  if (exon->objectType == CLASS_STICKYEXON) {
    stickyExon = (StickyExon *)exon;
    // sticky storing. Sticky exons contain normal exons ...

    for (i=0; i<StickyExon_getComponentExonCount(stickyExon); i++) {
      Exon *componentExon = StickyExon_getComponentExonAt(stickyExon,i);
      RawContig *contig; 

      if (!Exon_getContig(componentExon)->objectType != CLASS_RAWCONTIG) {
        fprintf(stderr,"Error: contig isn't raw contig when trying to store\n");
        exit(1);
      }
      contig = (RawContig *)Exon_getContig(componentExon);



      if (!contig || !RawContig_getDbID(contig)) {
        fprintf(stderr,"Component Exon does not have an attached contig "
                       "with a valid set database id. "
                       "Needs to have one set\n");
        exit(1);
      }

      if (!exonId) {
        sth->execute( sth2,
                      (IDType)RawContig_getDbID(contig),
                      Exon_getStart(componentExon),
                      Exon_getEnd(componentExon),
                      Exon_getStrand(componentExon),
                      Exon_getPhase(componentExon),
                      Exon_getEndPhase(componentExon),
                      Exon_getStickyRank(componentExon));
        exonId = sth->getInsertId(sth);
      } else {
        sth->execute( sth,
                      (IDType)exonId,
                      (IDType)RawContig_getDbID(contig),
                      Exon_getStart(componentExon),
                      Exon_getEnd(componentExon),
                      Exon_getStrand(componentExon),
                      Exon_getPhase(componentExon),
                      Exon_getEndPhase(componentExon),
                      Exon_getStickyRank(componentExon));
      }
    }
  } else {
    // normal storing
    RawContig *contig;

    if (!Exon_getContig(exon)->objectType != CLASS_RAWCONTIG) {
      fprintf(stderr,"Error: contig isn't raw contig when trying to store\n");
      exit(1);
    }
    contig = (RawContig *)Exon_getContig(exon);

    if (!contig || !RawContig_getDbID(contig)) {
      fprintf(stderr,"Exon does not have an attached contig with a valid " 
                     "database id.  Needs to have one set\n");
      exit(1);
    }

    sth->execute( sth2,
                  (IDType)RawContig_getDbID(contig),
                  Exon_getStart(exon),
                  Exon_getEnd(exon),
                  Exon_getStrand(exon),
                  Exon_getPhase(exon),
                  Exon_getEndPhase(exon),
                  Exon_getStickyRank(exon));
    exonId = sth->getInsertId(sth);
  }
  sth->finish(sth);

  if (Exon_getStableId(exon)) {
    if (!Exon_getCreated(exon) ||
        !Exon_getModified(exon) ||
        Exon_getVersion(exon) == -1) {
      fprintf(stderr, "Error: Trying to store incomplete stable id information for exon\n");
      exit(1);
    }

    sprintf(qStr,
        "INSERT INTO exon_stable_id(exon_id," 
        "version, stable_id, created, modified)"
                    " VALUES(" IDFMTSTR ",%d,'%s',FROM_UNIXTIME(%d),FROM_UNIXTIME(%d))",
         exonId,
         Exon_getVersion(exon),
         Exon_getStableId(exon),
         Exon_getCreated(exon),
         Exon_getModified(exon));


     sth = ea->prepare((BaseAdaptor *)ea,qStr,strlen(qStr));
     sth->execute(sth);
     sth->finish(sth);
   }


  // Now the supporting evidence
  // should be stored from featureAdaptor
  sprintf(qStr,
         "insert into supporting_feature (exon_id, feature_id, feature_type) "
         "values(%" IDFMTSTR ", %" IDFMTSTR ", %%s)");

  sth = ea->prepare((BaseAdaptor *)ea, qStr, strlen(qStr));

  dafa = DBAdaptor_getDNAAlignFeatureAdaptor(ea->dba);
  pafa = DBAdaptor_getProteinAlignFeatureAdaptor(ea->dba);

  if (stickyExon) {
    nExon = StickyExon_getComponentExonCount(stickyExon);
  } else {
    nExon = 1;
  }

/* nExon is 1 for non sticky and nComponent for sticky */
  for (i=0;i<nExon;i++) {
    int j;
    Exon *e;
    Vector *supportingFeatures;

    if (stickyExon) {
      e = StickyExon_getComponentExonAt(stickyExon,i);
    } else {
      e = exon;
    }

    supportingFeatures = Exon_getAllSupportingFeatures(e);

    for (j=0; j<Vector_getNumElement(supportingFeatures); j++) {
      BaseAlignFeature *sf = Vector_getElementAt(supportingFeatures, j);

      Class_assertType(CLASS_BASEALIGNFEATURE, sf->objectType);

      // sanity check
/* NIY
      if (!BaseAlignFeature_validate(sf)) {
        fprintf(stderr,"Warning: Supporting feature invalid. Skipping feature\n");
        continue;
      }
*/

      BaseAlignFeature_setContig(sf, Exon_getContig(e));

      if (Class_isDescendent(CLASS_DNADNAALIGNFEATURE, sf->objectType)) {
        DNAAlignFeatureAdaptor_store(dafa,sf);
        type = dnaType;
      } else if (Class_isDescendent(CLASS_DNAPEPALIGNFEATURE, sf->objectType)) {
        ProteinAlignFeatureAdaptor_store(pafa,sf);
        type = protType;
      } else {
        fprintf(stderr,"Warning: Supporting feature of unknown type. Skipping\n");
        continue;
      }

      sth->execute(sth, (IDType)exonId, BaseAlignFeature_getDbID(sf), type);
    }
  }
  sth->finish(sth);

  // 
  // Finally, update the dbID and adaptor of the exon (and any component exons)
  // to point to the new database
  // 

  if (stickyExon) {
    for (i=0; i<StickyExon_getComponentExonCount(stickyExon); i++) {
      Exon *e = StickyExon_getComponentExonAt(stickyExon,i);
      Exon_setDbID(e,exonId);
      Exon_setAdaptor(e,(BaseAdaptor *)ea);
    }
  }

  Exon_setAdaptor(exon,(BaseAdaptor *)ea);
  Exon_setDbID(exon, exonId);

  return Exon_getDbID(exon);
}

