#include "RepeatFeatureAdaptor.h"

#include "IDHash.h"
#include "RepeatConsensusAdaptor.h"
#include "AnalysisAdaptor.h"
#include "RawContigAdaptor.h"

#include "RepeatFeature.h"
#include "RepeatConsensus.h"

NameTableType RepeatFeatureAdaptor_tableNames = {{"repeat_feature","r"},
                                                 {"repeat_consensus","rc"},
                                                 {NULL,NULL}};

RepeatFeatureAdaptor *RepeatFeatureAdaptor_new(DBAdaptor *dba) {
  RepeatFeatureAdaptor *rfa;

  if ((rfa = (RepeatFeatureAdaptor *)calloc(1,sizeof(RepeatFeatureAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for RepeatFeatureAdaptor\n");
    return NULL;
  }
  BaseFeatureAdaptor_init((BaseFeatureAdaptor *)rfa, dba, REPEATFEATURE_ADAPTOR);

  rfa->getTables = RepeatFeatureAdaptor_getTables;
  rfa->getColumns = RepeatFeatureAdaptor_getColumns;
  rfa->store = RepeatFeatureAdaptor_store;
  rfa->objectsFromStatementHandle = RepeatFeatureAdaptor_objectsFromStatementHandle;
  rfa->defaultWhereClause = RepeatFeatureAdaptor_defaultWhereClause;

  return rfa;
}

int RepeatFeatureAdaptor_store(BaseFeatureAdaptor *bfa, Set *features) {
  RepeatConsensus *cons;
  RepeatConsensusAdaptor *rca = DBAdaptor_getRepeatConsensusAdaptor(bfa->dba);
  StatementHandle *sth;
  int i;
  char qStr[1024];

  sprintf(qStr,
    "INSERT into repeat_feature( repeat_feature_id"
                        ", contig_id"
                        ", contig_start"
                        ", contig_end"
                        ", contig_strand"
                        ", repeat_consensus_id"
                        ", repeat_start"
                        ", repeat_end"
                        ", score"
                        ", analysis_id )"
      " VALUES(NULL, %" IDFMTSTR ",%%d,%%d,%%d,%" IDFMTSTR ",%%d,%%d,%%f,%" IDFMTSTR ")");

  sth = bfa->prepare((BaseAdaptor *)bfa,qStr,strlen(qStr));

  for (i=0; i<Set_getNumElement(features); i++) {
    RepeatFeature *rf = Set_getElementAt(features,i);
    IDType dbID;
    IDType analId;
    IDType consId;
    RawContig *contig;
    RepeatConsensus *cons = RepeatFeature_getConsensus(rf);

    if (!cons) {
      fprintf(stderr,"Error: Must have a RepeatConsensus attached\n");
      exit(1);
    }
      
    // for tandem repeats - simply store consensus and repeat
    // one pair per hit. don't need to check consensi stored
    // already. consensus has name and class set to 'trf'

    if (!strcmp(RepeatConsensus_getRepeatClass(cons),"trf")) {

      // Look for matches already stored

// NIY This is a terribly slow way to do this - f**king consensi

      Set *match = RepeatConsensusAdaptor_fetchByClassAndSeq(rca, "trf", RepeatConsensus_getConsensus(cons)); 
      if (Set_getNumElement(match)) {
        RepeatConsensus *matchedCons = Set_getElementAt(match,0);
        RepeatConsensus_setDbID(cons,RepeatConsensus_getDbID(matchedCons));
        Set_free(match,RepeatConsensus_free);
      } else {
        RepeatConsensusAdaptor_store(rca,cons);
        Set_free(match,NULL);
      }

    } else if (!strcmp(RepeatConsensus_getRepeatClass(cons),"Simple_repeat")) {
      char tmpStr[EXTREMELEN];
      int len;
      // Look for matches already stored
      RepeatConsensus *match = RepeatConsensusAdaptor_fetchByNameAndClass(rca,RepeatConsensus_getName(cons),"Simple_repeat");

      strcpy(tmpStr,RepeatConsensus_getName(cons));
      len = strlen(tmpStr);

      if (tmpStr[0]=='(' &&
          tmpStr[len-2]==')' &&
          tmpStr[len-1]=='n') {
        memmove(tmpStr,&tmpStr[1],len);
        tmpStr[len-2] = '\0';
      }
      RepeatConsensus_setConsensus(cons,tmpStr);

      if (match) {
        RepeatConsensus_setDbID(cons,RepeatConsensus_getDbID(match));
        RepeatConsensus_free(match);
      } else {
        RepeatConsensusAdaptor_store(rca,cons);
      }
    } else {

      // for other repeats - need to see if a consensus is stored already
      if (!RepeatConsensus_getDbID(cons)) {
        RepeatConsensus *match = RepeatConsensusAdaptor_fetchByName(rca,RepeatConsensus_getName(cons));
    
        if (match) {
          //set the consensus dbID to be the same as the database one
          RepeatConsensus_setDbID(cons,RepeatConsensus_getDbID(match));
          RepeatConsensus_free(match);
        } else {
          // if we don't match a consensus already stored create a fake one 
          // and set consensus to 'N' as null seq not allowed
          // FIXME: not happy with this, but ho hum ...
          fprintf(stderr, "Warning: Can't find %s repeat consensus\n", RepeatConsensus_getName(cons));
          RepeatConsensus_setConsensus(cons,"N");
          RepeatConsensusAdaptor_store(rca,cons);
        }
      }
    }
    
    contig = RepeatFeature_getContig(rf);

/* NIY
    unless(defined $contig && $contig->isa("Bio::EnsEMBL::RawContig")) {
      $self->throw("RepeatFeature cannot be stored without a contig " .
           "attached via the attach_seq method");
    } unless($contig->dbID()) {
      $self->throw("RepeatFeature cannot be stored because attached contig " .
           "does not have a dbID");
    }
*/
    
    consId = RepeatConsensus_getDbID(RepeatFeature_getConsensus(rf));
    analId = Analysis_getDbID(RepeatFeature_getAnalysis(rf));
    sth->execute(sth,
          (IDType)RawContig_getDbID(contig),
          RepeatFeature_getStart(rf),
          RepeatFeature_getEnd(rf),
          RepeatFeature_getStrand(rf),
          (IDType)consId,
          RepeatFeature_getHitStart(rf),
          RepeatFeature_getHitEnd(rf),
          RepeatFeature_getScore(rf),
          (IDType)analId
         );

    dbID = sth->getInsertId(sth);

    RepeatFeature_setDbID(rf, dbID);
    RepeatFeature_setAdaptor(rf, (BaseAdaptor *)bfa);
  }
  sth->finish(sth);

  return 1;
}

NameTableType *RepeatFeatureAdaptor_getTables() {
  return &RepeatFeatureAdaptor_tableNames;
}

char *RepeatFeatureAdaptor_getColumns() {
  return "r.repeat_feature_id,"
         "r.contig_id,"
         "r.analysis_id,"
         "r.contig_start,"
         "r.contig_end,"
         "r.contig_strand,"
         "r.repeat_consensus_id,"
         "r.repeat_start,"
         "r.repeat_end,"
         "r.score,"
         "rc.repeat_name,"
         "rc.repeat_class,"
         "rc.repeat_consensus";
}

Set *RepeatFeatureAdaptor_objectsFromStatementHandle(BaseFeatureAdaptor *bfa,
                                                     StatementHandle *sth,
                                                     AssemblyMapper *mapper,
                                                     Slice *slice) {
  AnalysisAdaptor *aa;
  RawContigAdaptor *rca;
  RepeatConsensusAdaptor *rpca;
  Set *features;
  ResultRow *row;
  int i;
  IDHash *rcHash;

  aa = DBAdaptor_getAnalysisAdaptor(bfa->dba);
  rca = DBAdaptor_getRawContigAdaptor(bfa->dba);
  rpca = DBAdaptor_getRepeatConsensusAdaptor(bfa->dba);


  features = Set_new();
  rcHash = IDHash_new(IDHASH_SMALL);

  while (row = sth->fetchRow(sth)) {
    RepeatFeature *rf;
    Analysis  *analysis = AnalysisAdaptor_fetchByDbID(aa, row->getLongLongAt(row,2));
    RawContig *contig = RawContigAdaptor_fetchByDbID(rca, row->getLongLongAt(row,1));
    IDType repeatConsensusId = row->getLongLongAt(row,6);
    RepeatConsensus *rc;

    //create a repeat consensus object
    if (!IDHash_contains(rcHash, repeatConsensusId)) {
      rc = RepeatConsensus_new();
      RepeatConsensus_setDbID(rc, repeatConsensusId);
      RepeatConsensus_setName(rc, row->getStringAt(row,10));
      RepeatConsensus_setRepeatClass(rc, row->getStringAt(row,11));
      RepeatConsensus_setConsensus(rc, row->getStringAt(row,12));
      RepeatConsensus_setAdaptor(rc,(BaseAdaptor *)rpca);

      IDHash_add(rcHash,repeatConsensusId,rc);
    } else {
      rc = IDHash_getValue(rcHash,repeatConsensusId);
    }
    
    //create the new repeat feature
    rf = RepeatFeature_new();
    RepeatFeature_setDbID(rf,row->getLongLongAt(row,0));
    RepeatFeature_setStart(rf,row->getIntAt(row,3));
    RepeatFeature_setEnd(rf,row->getIntAt(row,4));
    RepeatFeature_setStrand(rf,row->getIntAt(row,5));
    RepeatFeature_setHitStart(rf,row->getIntAt(row,7));
    RepeatFeature_setHitEnd(rf,row->getIntAt(row,8));
    RepeatFeature_setContig(rf,contig); 
    RepeatFeature_setAnalysis(rf,analysis); 
    RepeatFeature_setConsensus(rf,rc); 

    if (row->col(row,9)) RepeatFeature_setScore(rf,row->getDoubleAt(row,9));

    Set_addElement(features,rf);
  }

  IDHash_free(rcHash,NULL);

  return features;
}

char *RepeatFeatureAdaptor_defaultWhereClause() {
  return "r.repeat_consensus_id = rc.repeat_consensus_id";
}

