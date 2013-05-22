/*
=head1 NAME

Bio::EnsEMBL::DBSQL::PredictionExonAdaptor - Performs database interaction for
PredictionExons.
*/
#include "PredictionExonAdaptor.h"

#include "PredictionTranscript.h"
#include "AnalysisAdaptor.h"
#include "SliceAdaptor.h"
#include "ChainedAssemblyMapper.h"
#include "AssemblyMapperAdaptor.h"



NameTableType PredictionExonAdaptor_tableNames = {{"prediction_exon","pe"},
                                                  {NULL,NULL}};

PredictionExonAdaptor *PredictionExonAdaptor_new(DBAdaptor *dba) {
  PredictionExonAdaptor *pea;

  if ((pea = (PredictionExonAdaptor *)calloc(1,sizeof(PredictionExonAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for PredictionExonAdaptor\n");
    return NULL;
  }
  BaseFeatureAdaptor_init((BaseFeatureAdaptor *)pea, dba, PREDICTIONEXON_ADAPTOR);

  pea->getTables                  = PredictionExonAdaptor_getTables;
  pea->getColumns                 = PredictionExonAdaptor_getColumns;
  pea->store                      = PredictionExonAdaptor_store;
  pea->objectsFromStatementHandle = PredictionExonAdaptor_objectsFromStatementHandle;
  pea->finalClause                = PredictionExonAdaptor_finalClause;

  return pea;
}


/*
#_tables
#
#  Arg [1]    : none
#  Example    : none
#  Description: PROTECTED implementation of superclass abstract method
#               returns the names, aliases of the tables to use for queries
#  Returntype : list of listrefs of strings
#  Exceptions : none
#  Caller     : internal
#
*/
NameTableType *PredictionExonAdaptor_getTables() {
  return &PredictionExonAdaptor_tableNames;
}

/*
#_columns
#
#  Arg [1]    : none
#  Example    : none
#  Description: PROTECTED implementation of superclass abstract method
#               returns a list of columns to use for queries
#  Returntype : list of strings
#  Exceptions : none
#  Caller     : internal
*/
char *PredictionExon_cols[] = {
                               "pe.prediction_exon_id",
                               "pe.seq_region_id",
                               "pe.seq_region_start",
                               "pe.seq_region_end",
                               "pe.seq_region_strand",
                               "pe.start_phase",
                               "pe.score",
                               "pe.p_value",
                               NULL
                              };

char **PredictionExonAdaptor_getColumns() {
  return PredictionExon_cols;
}


/*
# _final_clause
#
#  Arg [1]    : none
#  Example    : none
#  Description: PROTECTED implementation of superclass abstract method
#               returns a default end for the SQL-query ORDER BY)
#  Returntype : string
#  Exceptions : none
#  Caller     : internal
*/
char *PredictionExonAdaptor_finalClause() {
  return " ORDER BY pe.prediction_transcript_id, pe.exon_rank";
}


/*
=head2 fetch_all_by_PredictionTranscript

  Arg [1]    : Bio::EnsEMBL::PredcitionTranscript $transcript
  Example    : none
  Description: Retrieves all Exons for the Transcript in 5-3 order
  Returntype : listref Bio::EnsEMBL::Exon on Transcript slice 
  Exceptions : throws if transcript does not have a slice
  Caller     : Transcript->get_all_Exons()
  Status     : Stable

=cut
*/

Vector *PredictionExonAdaptor_fetchAllByPredictionTranscript(PredictionExonAdaptor *pea, PredictionTranscript *transcript) {
  // A mysterious option which, although referred to, doesn't seem to exist in the perl - 
  // use 'keep_all' option to keep exons that are off end of slice ???

  Slice *tSlice = PredictionTranscript_getSlice(transcript);
  Slice *slice;

  if (tSlice == NULL) {
    fprintf(stderr, "Transcript must have attached slice to retrieve exons.\n");
    exit(1);
  }

  // use a small slice the same size as the prediction transcript
  SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(pea->dba);
  slice = SliceAdaptor_fetchByFeature(sa, (SeqFeature *)transcript, 0, 0);

  char constraint[1024];
  sprintf(constraint, "pe.prediction_transcript_id = "IDFMTSTR, PredictionTranscript_getDbID(transcript));

  Vector *exons = PredictionExonAdaptor_fetchAllBySliceConstraint(pea, slice, constraint, NULL);

  // remap exon coordinates if necessary
  if (EcoString_strcmp(Slice_getName(slice), Slice_getName(tSlice))) {
    Vector *out = Vector_new();
    int i;
    for (i=0; i<Vector_getNumElement(exons); i++) {
      PredictionExon *ex = Vector_getElementAt(exons, i);
      Vector_addElement(out, PredictionExon_transfer(ex, tSlice));
    }
// NIY Do I need to set a free func here - probably???
    Vector_free(exons);

    exons = out;
  }

  return exons;
}



/*
=head2 store

  Arg [1]    : Bio::EnsEMBL::PredictionExon $exon
               The exon to store in this database
  Arg [2]    : int $prediction_transcript_id
               The internal identifier of the prediction exon that that this
               exon is associated with.
  Arg [3]    : int $rank
               The rank of the exon in the transcript (starting at 1)
  Example    : $pexon_adaptor->store($pexon, 1211, 2);
  Description: Stores a PredictionExon in the database
  Returntype : none
  Exceptions : thrown if exon does not have a slice attached
               or if $exon->start, $exon->end, $exon->strand, or $exon->phase 
               are not defined or if $exon is not a Bio::EnsEMBL::PredictionExon 
  Caller     : general
  Status     : Stable

=cut
*/
IDType PredictionExonAdaptor_store(PredictionExonAdaptor *pea, PredictionExon *pExon, IDType ptId, int rank) {
  fprintf(stderr,"PredictionExon store not implemented yet\n");
  exit(1);

/* NIY
  if(!ref($pexon) || !$pexon->isa('Bio::EnsEMBL::PredictionExon') ) {
    throw("Expected PredictionExon argument");
  }

  throw("Expected PredictionTranscript id argument.") if(!$pt_id);
  throw("Expected rank argument.") if(!$rank);

  my $db = $self->db();

  if($pexon->is_stored($db)) {
    warning('PredictionExon is already stored in this DB.');
    return $pexon->dbID();
  }

  if( ! $pexon->start || ! $pexon->end ||
      ! $pexon->strand || ! defined $pexon->phase ) {
    throw("PredictionExon does not have all attributes to store.\n" .
         "start, end, strand and phase attributes must be set.");
  }

  #maintain reference to original passed-in prediction exon
  my $original = $pexon;
  my $seq_region_id;
  ($pexon, $seq_region_id) = $self->_pre_store($pexon);

  my $sth = $db->dbc->prepare
    ("INSERT into prediction_exon (prediction_transcript_id, exon_rank, " .
                       "seq_region_id, seq_region_start, seq_region_end, " .
                       "seq_region_strand, start_phase, score, p_value) " .
      "VALUES ( ?, ?, ?, ?, ?, ?, ?, ?, ? )");

  $sth->bind_param(1,$pt_id,SQL_INTEGER);
  $sth->bind_param(2,$rank,SQL_SMALLINT);
  $sth->bind_param(3,$seq_region_id,SQL_INTEGER);
  $sth->bind_param(4,$pexon->start,SQL_INTEGER);
  $sth->bind_param(5,$pexon->end,SQL_INTEGER);
  $sth->bind_param(6,$pexon->strand,SQL_TINYINT);
  $sth->bind_param(7,$pexon->phase,SQL_TINYINT);
  $sth->bind_param(8,$pexon->score,SQL_DOUBLE);
  $sth->bind_param(9,$pexon->p_value,SQL_DOUBLE);

  $sth->execute();

  my $dbID = $sth->{'mysql_insertid'};

  #set the adaptor and dbID of the object they passed in
  $original->dbID($dbID);
  $original->adaptor($self);

  return $dbID;
*/
}



/* NIY
=head2 remove

  Arg [1]    : Bio::EnsEMBL::PredictionExon $exon
               the exon to remove from the database 
  Example    : $exon_adaptor->remove($exon);
  Description: Removes an exon from the database
  Returntype : none
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub remove {
  my $self = shift;
  my $pexon = shift;

  my $db = $self->db();

  if(!$pexon->is_stored($db)) {
    warning('PredictionExon is not in this DB - not removing');
    return undef;
  }

  my $sth = $self->prepare(
            "DELETE FROM prediction_exon WHERE prediction_exon_id = ?");
  $sth->bind_param( 1, $pexon->dbID, SQL_INTEGER );
  $sth->execute();

  $pexon->dbID(undef);
  $pexon->adaptor(undef);
}
*/



/*
=head2 list_dbIDs

  Arg [1]    : none
  Example    : @exon_ids = @{$exon_adaptor->list_dbIDs()};
  Description: Gets an array of internal ids for all exons in the current db
  Arg[1]     : <optional> int. not 0 for the ids to be sorted by the seq_region.
  Returntype : list of ints
  Exceptions : none
  Caller     : ?
  Status     : Stable

=cut
*/
Vector *PredictionExonAdaptor_listDbIDs(PredictionExonAdaptor *ea, int ordered) {

// NIY: Shouldn't really do this direct BaseAdaptor call but ...
  return BaseAdaptor_listDbIDs((BaseAdaptor *)ea, "prediction_exon", NULL, ordered);
}


/*
#_objs_from_sth

#  Arg [1]    : Hashreference $hashref
#  Example    : none 
#  Description: PROTECTED implementation of abstract superclass method.
#               responsible for the creation of Genes 
#  Returntype : listref of Bio::EnsEMBL::Genes in target coordinate system
#  Exceptions : none
#  Caller     : internal
#
*/
Vector *PredictionExonAdaptor_objectsFromStatementHandle(PredictionExonAdaptor *pea,
                                                         StatementHandle *sth,
                                                         AssemblyMapper *assMapper,
                                                         Slice *destSlice) {
  SliceAdaptor *sa     = DBAdaptor_getSliceAdaptor(pea->dba);

  Vector *exons = Vector_new();
  IDHash *sliceHash = IDHash_new(IDHASH_SMALL);

  long         destSliceStart;
  long         destSliceEnd;
  int          destSliceStrand;
  long         destSliceLength;
  CoordSystem *destSliceCs;
  char *       destSliceSrName;
  IDType       destSliceSrId = 0;
  AssemblyMapperAdaptor *asma;

  if (destSlice) {
    destSliceStart  = Slice_getStart(destSlice);
    destSliceEnd    = Slice_getEnd(destSlice);
    destSliceStrand = Slice_getStrand(destSlice);
    destSliceLength = Slice_getLength(destSlice);
    destSliceCs     = Slice_getCoordSystem(destSlice);
    destSliceSrName = Slice_getSeqRegionName(destSlice);
    destSliceSrId   = Slice_getSeqRegionId(destSlice);
    asma            = DBAdaptor_getAssemblyMapperAdaptor(pea->dba);
  }


  ResultRow *row;
  while (row = sth->fetchRow(sth)) {
    IDType predictionExonId = row->getLongLongAt(row, 0);
    IDType seqRegionId      = row->getLongLongAt(row, 1);
    long seqRegionStart     = row->getLongAt(row, 2);
    long seqRegionEnd       = row->getLongAt(row, 3);
    int seqRegionStrand     = row->getIntAt(row, 4);
    int startPhase          = row->getIntAt(row, 5);
    double score            = row->getDoubleAt(row, 6);
    double pValue           = row->getDoubleAt(row, 7);

    if (! IDHash_contains(sliceHash, seqRegionId)) {
      IDHash_add(sliceHash, seqRegionId, SliceAdaptor_fetchBySeqRegionId(sa, seqRegionId, POS_UNDEF, POS_UNDEF, STRAND_UNDEF));
    }
    Slice *slice = IDHash_getValue(sliceHash, seqRegionId);

    Slice *exonSlice = slice;

    char *srName      = Slice_getSeqRegionName(slice);
    CoordSystem *srCs = Slice_getCoordSystem(slice);


    // obtain a mapper if none was defined, but a dest_seq_region was
    if (assMapper == NULL &&
        destSlice != NULL &&
        CoordSystem_compare(destSliceCs, Slice_getCoordSystem(slice))) {
      assMapper = AssemblyMapperAdaptor_fetchByCoordSystems(asma, destSliceCs, Slice_getCoordSystem(slice));
    }
  
    //
    // Remap the feature coordinates to another coord system if a mapper
    // was provided.
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

      // Get a slice in the coord system we just mapped to
      //$slice = $slice_hash{ "ID:" . $seq_region_id } ||=
      //  $sa->fetch_by_seq_region_id($seq_region_id);
      if (! IDHash_contains(sliceHash, seqRegionId)) {
        IDHash_add(sliceHash, seqRegionId, SliceAdaptor_fetchBySeqRegionId(sa, seqRegionId, POS_UNDEF, POS_UNDEF, STRAND_UNDEF));
      }
      exonSlice = IDHash_getValue(sliceHash, seqRegionId);
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
      exonSlice = destSlice;
    }

    // Finally, create the new PredictionExon.
    PredictionExon *exon = PredictionExon_new();

    PredictionExon_setStart  (exon, seqRegionStart);
    PredictionExon_setEnd    (exon, seqRegionEnd);
    PredictionExon_setStrand (exon, seqRegionStrand);
    PredictionExon_setAdaptor(exon, (BaseAdaptor *)pea);
    PredictionExon_setSlice  (exon, exonSlice);
    PredictionExon_setDbID   (exon, predictionExonId);
    PredictionExon_setPhase  (exon, startPhase);
    PredictionExon_setScore  (exon, score);
    PredictionExon_setpValue (exon, pValue);

    Vector_addElement(exons, exon);
  }

  return exons;
}

