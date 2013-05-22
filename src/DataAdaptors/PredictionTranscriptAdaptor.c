/*
=head1 NAME

Bio::EnsEMBL::DBSQL::PredictionTranscriptAdaptor -
Performs database interaction related to PredictionTranscripts
*/

#include "PredictionTranscriptAdaptor.h"
#include "PredictionExonAdaptor.h"

#include "PredictionTranscript.h"
#include "PredictionExon.h"
#include "AnalysisAdaptor.h"
#include "SliceAdaptor.h"
#include "ChainedAssemblyMapper.h"

NameTableType PredictionTranscriptAdaptor_tableNames = {{"prediction_transcript","pt"},
                                                        {NULL,NULL}};

PredictionTranscriptAdaptor *PredictionTranscriptAdaptor_new(DBAdaptor *dba) {
  PredictionTranscriptAdaptor *pta;

  if ((pta = (PredictionTranscriptAdaptor *)calloc(1,sizeof(PredictionTranscriptAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for PredictionTranscriptAdaptor\n");
    return NULL;
  }
  BaseFeatureAdaptor_init((BaseFeatureAdaptor *)pta, dba, PREDICTIONTRANSCRIPT_ADAPTOR);

  pta->getTables = PredictionTranscriptAdaptor_getTables;
  pta->getColumns = PredictionTranscriptAdaptor_getColumns;
  pta->store = PredictionTranscriptAdaptor_store;
  pta->objectsFromStatementHandle = PredictionTranscriptAdaptor_objectsFromStatementHandle;

  return pta;
}

/*
# _tables
#
#  Arg [1]    : none
#  Example    : none
#  Description: Implements abstract superclass method to define the table used
#               to retrieve prediction transcripts from the database
#  Returntype : string
#  Exceptions : none
#  Caller     : generic_fetch
*/
NameTableType *PredictionTranscriptAdaptor_getTables(void) {
  return &PredictionTranscriptAdaptor_tableNames;
}


/*
# _columns

#  Arg [1]    : none
#  Example    : none
#  Description: Implements abstract superclass method to define the columns
#               retrieved in database queries used to create prediction 
#               transcripts.
#  Returntype : list of strings
#  Exceptions : none
#  Caller     : generic_fetch
#
*/

char *PredictionTranscript_cols[] = { 
                                     "pt.prediction_transcript_id",
                                     "pt.seq_region_id",
                                     "pt.seq_region_start",
                                     "pt.seq_region_end",
                                     "pt.seq_region_strand",
                                     "pt.analysis_id",
                                     "pt.display_label",
                                     NULL };
char **PredictionTranscriptAdaptor_getColumns(void) {
  return PredictionTranscript_cols;
}


/*
=head2 fetch_by_stable_id

  Arg [1]    : string $stable_id
               The stable id of the transcript to retrieve
  Example    : $trans = $trans_adptr->fetch_by_stable_id('GENSCAN00000001234');
  Description: Retrieves a prediction transcript via its display_label.
               This method is called fetch_by_stable_id for polymorphism with
               the TranscriptAdaptor.  Prediction transcript display_labels are
               not necessarily stable in that the same identifier may be reused
               for a completely different prediction transcript in a subsequent
               database release.
  Returntype : Bio::EnsEMBL::PredictionTranscript
  Caller     : general
  Status     : Stable

=cut
*/

PredictionTranscript *PredictionTranscriptAdaptor_fetchByStableId(PredictionTranscriptAdaptor *pta, char *stableId) {
  if (stableId == NULL) {
    fprintf(stderr,"Error: Stable_id argument expected in PredictionTranscriptAdaptor_fetchByStableId\n");
    exit(1);
  }

  NameTableType *tables = pta->getTables();
  char **primTab = (*tables)[0];
  char *tableSynonym = primTab[SYN];

  char constraint[1024];
  sprintf(constraint, "%s.display_label = '%s'", tableSynonym, stableId);

  Vector *pts = PredictionTranscriptAdaptor_genericFetch(pta, constraint, NULL, NULL);
  PredictionTranscript *pt = Vector_getElementAt(pts, 0);
  Vector_free(pts);
// NIY: Free pts if there are more than 1

// NIY: Perl seemed to allow no pts so undef return, but I've not done that in other places eg. TranscriptAdaptor so didn't here either
  return pt;
}

/*

=head2 fetch_all_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               The slice to fetch transcripts on.
  Arg [3]    : (optional) boolean $load_exons
               if true, exons will be loaded immediately rather than
               lazy loaded later.
  Example    : $transcripts = $
  Description: Overrides superclass method to optionally load exons
               immediately rather than lazy-loading them later.  This
               is more efficient when there are a lot of transcripts whose
               exons are going to be used.
  Returntype : reference to list of transcripts
  Exceptions : thrown if exon cannot be placed on transcript slice
  Caller     : Slice::get_all_Transcripts
  Status     : Stable

=cut
*/

typedef struct PredictionTranscriptRankPairStruct {
  PredictionTranscript *transcript;
  int rank;
} PredictionTranscriptRankPair;

PredictionTranscriptRankPair *PredictionTranscriptRankPair_new(PredictionTranscript *transcript, int rank) {
  PredictionTranscriptRankPair *trp;

  if ((trp = (PredictionTranscriptRankPair *)calloc(1,sizeof(PredictionTranscriptRankPair))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for PredictionTranscriptRankPair\n");
    return NULL;
  }
  trp->transcript = transcript;
  trp->rank = rank;

  return trp;
}

void PredictionTranscriptRankPair_free(PredictionTranscriptRankPair *trp) {
  free(trp);
}

// This is almost entirely a copy and paste from Transcript

Vector *PredictionTranscriptAdaptor_fetchAllBySlice(PredictionTranscriptAdaptor *pta, Slice *slice, char *logicName, int loadExons) {

  //my $transcripts = $self->SUPER::fetch_all_by_Slice($slice,$logic_name);
  Vector *transcripts = BaseFeatureAdaptor_fetchAllBySlice((BaseFeatureAdaptor *)pta, slice, logicName);

  // if there are 0 or 1 transcripts still do lazy-loading
  if ( ! loadExons || Vector_getNumElement(transcripts) < 2 ) {
    return transcripts;
  }

  // preload all of the exons now, instead of lazy loading later
  // faster than 1 query per transcript

  // get extent of region spanned by transcripts
  long minStart =  2000000000;
  long maxEnd   = -2000000000;

  int i;
  for (i=0; i<Vector_getNumElement(transcripts); i++) {
    PredictionTranscript *t  = Vector_getElementAt(transcripts, i);
    if (PredictionTranscript_getSeqRegionStart(t) < minStart) {
      minStart = PredictionTranscript_getSeqRegionStart(t);
    }
    if (PredictionTranscript_getSeqRegionEnd(t) > maxEnd) {
      maxEnd = PredictionTranscript_getSeqRegionEnd(t);
    }
  }

  Slice *extSlice;

  if (minStart >= Slice_getStart(slice) && maxEnd <= Slice_getEnd(slice)) {
    extSlice = slice;
  } else {
    SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(pta->dba);
    extSlice = SliceAdaptor_fetchByRegion(sa, Slice_getCoordSystemName(slice), Slice_getSeqRegionName(slice),
                                          minStart, maxEnd, Slice_getStrand(slice), CoordSystem_getVersion(Slice_getCoordSystem(slice)), 0);
  }

  // associate exon identifiers with transcripts
  IDHash *trHash = IDHash_new(IDHASH_MEDIUM);
  for (i=0; i<Vector_getNumElement(transcripts); i++) {
    PredictionTranscript *t  = Vector_getElementAt(transcripts, i);
    if ( ! IDHash_contains(trHash, PredictionTranscript_getDbID(t))) {
      IDHash_add(trHash, PredictionTranscript_getDbID(t), t);
    }
  }

  IDType *uniqueIds = IDHash_getKeys(trHash);

  char tmpStr[1024];
  char qStr[655500];
  int lenNum;
  int endPoint = sprintf(qStr, "SELECT prediction_transcript_id, prediction_exon_id, exon_rank FROM prediction_exon WHERE  prediction_transcript_id IN (");
  for (i=0; i<IDHash_getNumValues(trHash); i++) {
    if (i!=0) {
      qStr[endPoint++] = ',';
      qStr[endPoint++] = ' ';
    }
    lenNum = sprintf(tmpStr,IDFMTSTR,uniqueIds[i]);
    memcpy(&(qStr[endPoint]), tmpStr, lenNum);
    endPoint+=lenNum;
  }
  qStr[endPoint++] = ')';
  qStr[endPoint] = '\0';

  free(uniqueIds);

  StatementHandle *sth = pta->prepare((BaseAdaptor *)pta,qStr,strlen(qStr));
  sth->execute(sth);

  IDHash *exTrHash = IDHash_new(IDHASH_MEDIUM);
  ResultRow *row;
  while (row = sth->fetchRow(sth)) {
    IDType trId = row->getLongLongAt(row,0);
    IDType exId = row->getLongLongAt(row,1);
    int    rank = row->getIntAt(row,2);

    if (! IDHash_contains(exTrHash, exId)) {
      Vector *vec = Vector_new();
      Vector_setFreeFunc(vec, PredictionTranscriptRankPair_free);
      IDHash_add(exTrHash, exId, vec);
    }
    Vector *exVec = IDHash_getValue(exTrHash, exId);
    PredictionTranscriptRankPair *trp = PredictionTranscriptRankPair_new(IDHash_getValue(trHash, trId), rank);
    Vector_addElement(exVec, trp);
  }

  IDHash_free(trHash, NULL);

  sth->finish(sth);

  PredictionExonAdaptor *pea = DBAdaptor_getPredictionExonAdaptor(pta->dba);
  Vector *exons = PredictionExonAdaptor_fetchAllBySlice(pea, extSlice);

  // move exons onto transcript slice, and add them to transcripts
  for (i=0; i<Vector_getNumElement(exons); i++) {
    PredictionExon *ex = Vector_getElementAt(exons, i);

  // Perl didn't have this line - it was in GeneAdaptor version so I think I'm going to keep it
    if (!IDHash_contains(exTrHash, PredictionExon_getDbID(ex))) continue;

    Exon *newEx;
    if (slice != extSlice) {
      newEx = PredictionExon_transfer(ex, slice);
      if (newEx == NULL) {
        fprintf(stderr, "Unexpected. Exon could not be transferred onto PredictionTranscript slice.\n");
        exit(1);
      }
    } else {
      newEx = ex;
    }

    Vector *exVec = IDHash_getValue(exTrHash, PredictionExon_getDbID(newEx));
    int j;
    for (j=0; j<Vector_getNumElement(exVec); j++) {
      PredictionTranscriptRankPair *trp = Vector_getElementAt(exVec, j);
      PredictionTranscript_addExon(trp->transcript, newEx, &trp->rank);
    }
  }

  IDHash_free(exTrHash, Vector_free);

  return transcripts;
}


/*
=head2 _objs_from_sth

  Arg [1]    : DBI:st $sth 
               An executed DBI statement handle
  Arg [2]    : (optional) Bio::EnsEMBL::Mapper $mapper 
               An mapper to be used to convert contig coordinates
               to assembly coordinates.
  Arg [3]    : (optional) Bio::EnsEMBL::Slice $slice
               A slice to map the prediction transcript to.   
  Example    : $p_transcripts = $self->_objs_from_sth($sth);
  Description: Creates a list of Prediction transcripts from an executed DBI
               statement handle.  The columns retrieved via the statement 
               handle must be in the same order as the columns defined by the
               _columns method.  If the slice argument is provided then the
               the prediction transcripts will be in returned in the coordinate
               system of the $slice argument.  Otherwise the prediction 
               transcripts will be returned in the RawContig coordinate system.
  Returntype : reference to a list of Bio::EnsEMBL::PredictionTranscripts
  Exceptions : none
  Caller     : superclass generic_fetch
  Status     : Stable

=cut
*/
Vector *PredictionTranscriptAdaptor_objectsFromStatementHandle(PredictionTranscriptAdaptor *pta, 
                                                               StatementHandle *sth, 
                                                               AssemblyMapper *assMapper, 
                                                               Slice *destSlice) {
  SliceAdaptor *sa     = DBAdaptor_getSliceAdaptor(pta->dba);
  AnalysisAdaptor *aa  = DBAdaptor_getAnalysisAdaptor(pta->dba);

  Vector *pTranscripts = Vector_new();
  IDHash *sliceHash = IDHash_new(IDHASH_SMALL);

  long         destSliceStart;
  long         destSliceEnd;
  int          destSliceStrand;
  long         destSliceLength;
  char *       destSliceSrName;
  IDType       destSliceSrId = 0;

  if (destSlice) {
    destSliceStart  = Slice_getStart(destSlice);
    destSliceEnd    = Slice_getEnd(destSlice);
    destSliceStrand = Slice_getStrand(destSlice);
    destSliceLength = Slice_getLength(destSlice);
    destSliceSrName = Slice_getSeqRegionName(destSlice);
    destSliceSrId   = Slice_getSeqRegionId(destSlice);
  }

  ResultRow *row;
  while (row = sth->fetchRow(sth)) {
    IDType predictionTranscriptId = row->getLongLongAt(row,0);
    IDType seqRegionId            = row->getLongLongAt(row,1);
    long seqRegionStart           = row->getLongAt(row,2);
    long seqRegionEnd             = row->getLongAt(row,3);
    int seqRegionStrand           = row->getIntAt(row,4);
    IDType analysisId             = row->getLongLongAt(row,5);
    char *displayLabel            = row->getStringAt(row,6);

    // get the analysis object
    Analysis *analysis = AnalysisAdaptor_fetchByDbID(aa, analysisId);

    if (! IDHash_contains(sliceHash, seqRegionId)) {
      IDHash_add(sliceHash, seqRegionId, SliceAdaptor_fetchBySeqRegionId(sa, seqRegionId, POS_UNDEF, POS_UNDEF, STRAND_UNDEF));
    }
    Slice *slice = IDHash_getValue(sliceHash, seqRegionId);

    Slice *ptSlice = slice;

    char *srName      = Slice_getSeqRegionName(slice);
    CoordSystem *srCs = Slice_getCoordSystem(slice);

    //
    // remap the feature coordinates to another coord system
    // if a mapper was provided
    //
    if (assMapper != NULL) {
      MapperRangeSet *mrs;

      // Slightly suspicious about need for this if statement so left in perl statements for now
      if (destSlice != NULL &&
          assMapper->objectType == CLASS_CHAINEDASSEMBLYMAPPER) {
        MapperRangeSet *mrs = ChainedAssemblyMapper_map(assMapper, srName, seqRegionStart, seqRegionEnd, seqRegionStrand, srCs, 1, destSlice);
      } else {
        MapperRangeSet *mrs = AssemblyMapper_fastMap(assMapper, srName, seqRegionStart, seqRegionEnd, seqRegionStrand, srCs, NULL);
      }

      // skip features that map to gaps or coord system boundaries
      if (MapperRangeSet_getNumRange(mrs) == 0) {
        continue;
      }
      MapperRange *range = MapperRangeSet_getRangeAt(mrs, 0);
      if (range->rangeType == MAPPERRANGE_GAP) {
        fprintf(stderr,"Got a mapper gap in gene obj_from_sth - not sure if this is allowed\n");
        exit(1);
      } else {
        MapperCoordinate *mc = (MapperCoordinate *)range;

        seqRegionId     = mc->id;
        seqRegionStart  = mc->start;
        seqRegionEnd    = mc->end;
        seqRegionStrand = mc->strand;
      }

      MapperRangeSet_free(mrs);

      if (! IDHash_contains(sliceHash, seqRegionId)) {
        IDHash_add(sliceHash, seqRegionId, SliceAdaptor_fetchBySeqRegionId(sa, seqRegionId, POS_UNDEF, POS_UNDEF, STRAND_UNDEF));
      }
      ptSlice = IDHash_getValue(sliceHash, seqRegionId);
    }

    //
    // If a destination slice was provided convert the coords
    // If the dest_slice starts at 1 and is foward strand, nothing needs doing
    //
    if (destSlice != NULL) {
      if (destSliceStart != 1 || destSliceStrand != 1) {
        if (destSliceStrand == 1) {
          seqRegionStart = seqRegionStart - destSliceStart + 1;
          seqRegionEnd   = seqRegionEnd - destSliceStart + 1;
        } else {
          long tmpSeqRegionStart = seqRegionStart;
          seqRegionStart = destSliceEnd - seqRegionEnd + 1;
          seqRegionEnd   = destSliceEnd - tmpSeqRegionStart + 1;

          seqRegionStrand = -seqRegionStrand;
        }
      }
      // throw away features off the end of the requested slice
      if (seqRegionEnd < 1 || seqRegionStart > destSliceLength || (destSliceSrId != seqRegionId)) {
        continue;
      }
      ptSlice = destSlice;
    }
    
    // Finally, create the new PredictionTranscript.
    PredictionTranscript *pt = PredictionTranscript_new();

    PredictionTranscript_setStart       (pt, seqRegionStart);
    PredictionTranscript_setEnd         (pt, seqRegionEnd);
    PredictionTranscript_setStrand      (pt, seqRegionStrand);
    PredictionTranscript_setSlice       (pt, ptSlice);
    PredictionTranscript_setAnalysis    (pt, analysis);
    PredictionTranscript_setAdaptor     (pt, (BaseAdaptor *)pta);
    PredictionTranscript_setDbID        (pt, predictionTranscriptId);
    PredictionTranscript_setDisplayLabel(pt, displayLabel);

    Vector_addElement(pTranscripts, pt);
  }

  return pTranscripts;
}


/* NIY
=head2 store

  Arg [1]    : list of Bio::EnsEMBL::PredictionTranscript @pre_transcripts 
  Example    : $prediction_transcript_adaptor->store(@pre_transcripts);
  Description: Stores a list of given prediction transcripts in database. 
               Puts dbID and Adaptor into each object stored object.
  Returntype : none
  Exceptions : on wrong argument type 
  Caller     : general 
  Status     : Stable

=cut
*/

void PredictionTranscriptAdaptor_store(PredictionTranscriptAdaptor *pta, Vector *preTranscripts) {
  fprintf(stderr,"PredictionTranscript store not implemented yet\n");
  exit(1);


/*
  my $ptstore_sth = $self->prepare
    (qq{INSERT INTO prediction_transcript (seq_region_id, seq_region_start,
                                           seq_region_end, seq_region_strand, 
                                           analysis_id, display_label)
        VALUES( ?, ?, ?, ?, ?, ?)});

  my $ptupdate_sth = $self->prepare
    (qq{UPDATE prediction_transcript SET display_label = ?
        WHERE  prediction_transcript_id = ?});

  my $db = $self->db();
  my $analysis_adaptor = $db->get_AnalysisAdaptor();
  my $pexon_adaptor = $db->get_PredictionExonAdaptor();

  FEATURE: foreach my $pt (@pre_transcripts) {
    if(!ref($pt) || !$pt->isa('Bio::EnsEMBL::PredictionTranscript')) {
      throw('Expected PredictionTranscript argument not [' . ref($pt).']');
    }

    #skip prediction transcripts that have already been stored
    if($pt->is_stored($db)) {
      warning('Not storing already stored prediction transcript '. $pt->dbID);
      next FEATURE;
    }

    #get analysis and store it if it is not in the db
    my $analysis = $pt->analysis();
    if(!$analysis) {
      throw('Prediction transcript must have analysis to be stored.');
    }
    if(!$analysis->is_stored($db)) {
      $analysis_adaptor->store($analysis);
    }

    #ensure that the transcript coordinates are correct, they may not be,
    #if somebody has done some exon coordinate juggling and not recalculated
    #the transcript coords.
    $pt->recalculate_coordinates();

    my $original = $pt;
    my $seq_region_id;
    ($pt, $seq_region_id) = $self->_pre_store($pt);

    #store the prediction transcript
    $ptstore_sth->bind_param(1,$seq_region_id,SQL_INTEGER);
    $ptstore_sth->bind_param(2,$pt->start,SQL_INTEGER);
    $ptstore_sth->bind_param(3,$pt->end,SQL_INTEGER);
    $ptstore_sth->bind_param(4,$pt->strand,SQL_TINYINT);
    $ptstore_sth->bind_param(5,$analysis->dbID,SQL_INTEGER);
    $ptstore_sth->bind_param(6,$pt->display_label,SQL_VARCHAR);

    $ptstore_sth->execute();

    my $pt_id = $ptstore_sth->{'mysql_insertid'};
    $original->dbID($pt_id);
    $original->adaptor($self);

    #store the exons
    my $rank = 1;
    foreach my $pexon (@{$original->get_all_Exons}) {
      $pexon_adaptor->store($pexon, $pt_id, $rank++);
    }

    # if a display label was not defined autogenerate one
    if(!defined($pt->display_label())) {
      my $zeros = '0' x (11 - length($pt_id));
      my $display_label = uc($analysis->logic_name()) . $zeros . $pt_id;
      $ptupdate_sth->bind_param(1,$display_label,SQL_VARCHAR);
      $ptupdate_sth->bind_param(2,$pt_id,SQL_INTEGER);
      $ptupdate_sth->execute();
      $original->display_label($display_label);
    }
  }
*/
}



/* NIY
=head2 remove

  Arg [1]    : Bio::EnsEMBL::PredictionTranscript $pt 
  Example    : $prediction_transcript_adaptor->remove($pt);
  Description: removes given prediction transcript $pt from database. 
  Returntype : none
  Exceptions : throws if argument not a  Bio::EnsEMBL::PredictionTranscript
  Caller     : general
  Status     : Stable

=cut

sub remove {
  my $self = shift;
  my $pre_trans = shift;

  if(!ref($pre_trans)||!$pre_trans->isa('Bio::EnsEMBL::PredictionTranscript')){
    throw('Expected PredictionTranscript argument.');
  }

  if(!$pre_trans->is_stored($self->db())) {
    warning('PredictionTranscript is not stored in this DB - not removing.');
    return;
  }

  #remove all associated prediction exons
  my $pexon_adaptor = $self->db()->get_PredictionExonAdaptor();
  foreach my $pexon (@{$pre_trans->get_all_Exons}) {
    $pexon_adaptor->remove($pexon);
  }

  #remove the prediction transcript
  my $sth = $self->prepare( "DELETE FROM prediction_transcript
                             WHERE prediction_transcript_id = ?" );
  $sth->bind_param(1,$pre_trans->dbID,SQL_INTEGER);
  $sth->execute();

  #unset the adaptor and internal id
  $pre_trans->dbID(undef);
  $pre_trans->adaptor(undef);
}
*/


/*
=head2 list_dbIDs

  Arg [1]    : none
  Example    : @feature_ids = @{$prediction_transcript_adaptor->list_dbIDs()};
  Description: Gets an array of internal ids for all prediction transcript
               features in the current db
  Arg[1]     : <optional> int. not 0 for the ids to be sorted by the seq_region.
  Returntype : list of ints
  Exceptions : none
  Caller     : ?
  Status     : Stable

=cut
*/
Vector *PredictionTranscriptAdaptor_listDbIDs(PredictionTranscriptAdaptor *pta, int ordered) {

// NIY: Shouldn't really do this direct BaseAdaptor call but ...
  return BaseAdaptor_listDbIDs((BaseAdaptor *)pta, "prediction_transcript", NULL, ordered);
}
