#include "RepeatFeatureAdaptor.h"

#include "IDHash.h"
#include "RepeatConsensusAdaptor.h"
#include "AnalysisAdaptor.h"
#include "SliceAdaptor.h"

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

int RepeatFeatureAdaptor_store(BaseFeatureAdaptor *bfa, Vector *features) {
  RepeatConsensusAdaptor *rca = DBAdaptor_getRepeatConsensusAdaptor(bfa->dba);
  StatementHandle *sth;
  int i;
  char qStr[1024];

  sprintf(qStr,
    "INSERT into repeat_feature( repeat_feature_id"
                        ", seq_region_id"
                        ", seq_region_start"
                        ", seq_region_end"
                        ", seq_region_strand"
                        ", repeat_consensus_id"
                        ", repeat_start"
                        ", repeat_end"
                        ", score"
                        ", analysis_id )"
      " VALUES(NULL, %" IDFMTSTR ",%%d,%%d,%%d,%" IDFMTSTR ",%%d,%%d,%%f,%" IDFMTSTR ")");

  sth = bfa->prepare((BaseAdaptor *)bfa,qStr,strlen(qStr));

  for (i=0; i<Vector_getNumElement(features); i++) {
    RepeatFeature *rf = Vector_getElementAt(features,i);
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

      Vector *match = RepeatConsensusAdaptor_fetchByClassAndSeq(rca, "trf", RepeatConsensus_getConsensus(cons)); 
      Vector_setFreeFunc(match,Object_freeImpl);
      if (Vector_getNumElement(match)) {
        RepeatConsensus *matchedCons = Vector_getElementAt(match,0);
        RepeatConsensus_setDbID(cons,RepeatConsensus_getDbID(matchedCons));
        Vector_free(match);
      } else {
        Vector *consVector = Vector_new();
        Vector_addElement(consVector,cons);
        RepeatConsensusAdaptor_store(rca,consVector);
        Vector_free(match);
        Vector_free(consVector);
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
        Vector *consVector = Vector_new();
        Vector_addElement(consVector,cons);
        RepeatConsensusAdaptor_store(rca,consVector);
        Vector_free(consVector);
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
          Vector *consVector = Vector_new();
          // if we don't match a consensus already stored create a fake one 
          // and set consensus to 'N' as null seq not allowed
          // FIXME: not happy with this, but ho hum ...
          fprintf(stderr, "Warning: Can't find %s repeat consensus\n", RepeatConsensus_getName(cons));
          RepeatConsensus_setConsensus(cons,"N");
          Vector_addElement(consVector,cons);
          RepeatConsensusAdaptor_store(rca,consVector);
          Vector_free(consVector);
        }
      }
    }
    
    if (RepeatFeature_getContig(rf)->objectType != CLASS_RAWCONTIG) {
      fprintf(stderr,"Error: contig isn't raw contig when trying to store (is %d)\n",RepeatFeature_getContig(rf)->objectType);
      exit(1);
    }

    contig = (RawContig *)RepeatFeature_getContig(rf);

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

char *RepeatFeature_cols[] = {
         "r.repeat_feature_id",
         "r.seq_region_id",
         "r.analysis_id",
         "r.seq_region_start",
         "r.seq_region_end",
         "r.seq_region_strand",
         "r.repeat_consensus_id",
         "r.repeat_start",
         "r.repeat_end",
         "r.score",
         "rc.repeat_name",
         "rc.repeat_class",
         "rc.repeat_consensus",
         NULL };

char **RepeatFeatureAdaptor_getColumns() {
  return RepeatFeature_cols;
}

/*
=head2 fetch_all_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
  Arg [2]    : (optional) string $logic_name
               Limits RepeatFeatures obtained to those having an Analysis with
               of the specified logic_name.  If no logic name is specified
               Repeats of all analysis types are retrieved.
  Arg [3]    : (optional) string/array $repeat_type
               Limits RepeatFeatures obtained to those of specified
               repeat_type
  Example    : @rfeats = @{$rfa->fetch_all_by_Slice($slice, undef, 'Type II Transposons')};
               @rfeats = @{$rfa->fetch_all_by_Slice($slice, undef, ['Type II Transposons', 'RNA repeats'])};
  Description: Retrieves repeat features overlapping the area designated by
               the provided slice argument.  Returned features will be in
               in the same coordinate system as the provided slice and will
               have coordinates relative to the slice start.
  Returntype : reference to a list of Bio::EnsEMBL::RepeatFeatures.
  Exceptions : throw on bad argument
  Caller     : Slice::get_all_RepeatFeatures
  Status     : Stable

=cut
*/
// New!!!
Vector *RepeatFeatureAdaptor_fetchAllBySlice(RepeatFeatureAdaptor *rfa, Slice *slice, char *logicName, Vector *repeatTypes) {
  char constraint[1024];
  constraint[0] = '\0';

  // MySQL was optimising the query the incorrect way when joining to
  // the repeat_consensus table on type

  // Hack - direct access
  rfa->straightJoinFlag = 1;

  if (repeatTypes != NULL && Vector_getNumElement(repeatTypes) > 0) {
    if (Vector_getNumElement(repeatTypes) > 1) {
      strcpy(constraint, "rc.repeat_type IN (");
      int i;
      for (i=0; i<Vector_getNumElement(repeatTypes); i++) {
        char *repeatType = Vector_getElementAt(repeatTypes, i);
        if (i!=0) {
          strcat(constraint, ", ");
        }
        sprintf(constraint, "'%s'", repeatType);
      }
      strcat(constraint,")");
    } else {
      char *repeatType = Vector_getElementAt(repeatTypes, 0);
      sprintf(constraint, "rc.repeat_type = '%s'", repeatType);
    }
  }

  Vector *result = RepeatFeatureAdaptor_fetchAllBySliceConstraint(rfa, slice, constraint, logicName);

  rfa->straightJoinFlag = 0;

  return result;
}

Vector *RepeatFeatureAdaptor_objectsFromStatementHandle(BaseFeatureAdaptor *bfa,
                                                     StatementHandle *sth,
                                                     AssemblyMapper *mapper,
                                                     Slice *destSlice) {
  AnalysisAdaptor *aa;
  SliceAdaptor *sa;
  RepeatConsensusAdaptor *rpca;
  Vector *features;
  ResultRow *row;
  IDHash *rcHash;

  aa = DBAdaptor_getAnalysisAdaptor(bfa->dba);
  sa = DBAdaptor_getSliceAdaptor(bfa->dba);
  rpca = DBAdaptor_getRepeatConsensusAdaptor(bfa->dba);


  features = Vector_new();
  rcHash = IDHash_new(IDHASH_SMALL);

  while ((row = sth->fetchRow(sth))) {
    RepeatFeature *rf;
    Analysis  *analysis = AnalysisAdaptor_fetchByDbID(aa, row->getLongLongAt(row,2));
    Slice *slice = SliceAdaptor_fetchBySeqRegionId(sa, row->getLongLongAt(row,1), POS_UNDEF, POS_UNDEF, 1);
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
    RepeatFeature_setContig(rf,slice); 
    RepeatFeature_setAnalysis(rf,analysis); 
    RepeatFeature_setConsensus(rf,rc); 

    if (row->col(row,9)) RepeatFeature_setScore(rf,row->getDoubleAt(row,9));

    Vector_addElement(features,rf);
  }

  IDHash_free(rcHash,NULL);

  return features;
}

char *RepeatFeatureAdaptor_defaultWhereClause() {
  return "r.repeat_consensus_id = rc.repeat_consensus_id";
}

