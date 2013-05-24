/*
=head1 NAME

Bio::EnsEMBL::DBSQL::SimpleFeatureAdaptor
*/

#include "SimpleFeatureAdaptor.h"

#include "SliceAdaptor.h"
#include "AnalysisAdaptor.h"
#include "ChainedAssemblyMapper.h"

NameTableType SimpleFeatureAdaptor_tableNames = {{"simple_feature","sf"},{NULL,NULL}};

SimpleFeatureAdaptor *SimpleFeatureAdaptor_new(DBAdaptor *dba) {
  SimpleFeatureAdaptor *sfa;

  if ((sfa = (SimpleFeatureAdaptor *)calloc(1,sizeof(SimpleFeatureAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for SimpleFeatureAdaptor\n");
    return NULL;
  }
  BaseFeatureAdaptor_init((BaseFeatureAdaptor *)sfa, dba, SIMPLEFEATURE_ADAPTOR);

  sfa->getTables = SimpleFeatureAdaptor_getTables;
  sfa->getColumns = SimpleFeatureAdaptor_getColumns;
  sfa->store = SimpleFeatureAdaptor_store;
  sfa->objectsFromStatementHandle = SimpleFeatureAdaptor_objectsFromStatementHandle;


  return sfa;
}


/*
=head2 _tables

  Arg [1]    : none
  Example    : none
  Description: PROTECTED implementation of superclass abstract method
               returns the names, aliases of the tables to use for queries
  Returntype : list of listrefs of strings
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut
*/
NameTableType *SimpleFeatureAdaptor_getTables(void) {
  return &SimpleFeatureAdaptor_tableNames;
}

/*
=head2 _columns

  Arg [1]    : none
  Example    : none
  Description: PROTECTED implementation of superclass abstract method
               returns a list of columns to use for queries
  Returntype : list of strings
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut
*/

char *SimpleFeature_cols[] =
        {"sf.simple_feature_id",
         "sf.seq_region_id",
         "sf.seq_region_start",
         "sf.seq_region_end",
         "sf.seq_region_strand",
         "sf.display_label",
         "sf.analysis_id",
         "sf.score",
         NULL};

char **SimpleFeatureAdaptor_getColumns(void) {
  return SimpleFeature_cols;
}


/*
=head2 store

  Arg [1]    : list of Bio::EnsEMBL::SimpleFeatures @sf
               the simple features to store in the database
  Example    : $simple_feature_adaptor->store(@simple_feats);
  Description: Stores a list of simple feature objects in the database
  Returntype : none
  Exceptions : thrown if @sf is not defined, if any of the features do not
               have an attached slice.
               or if any elements of @sf are not Bio::EnsEMBL::SimpleFeatures 
  Caller     : general
  Status     : Stable

=cut
*/

int SimpleFeatureAdaptor_store(BaseFeatureAdaptor *bfa, Vector *features) {
  fprintf(stderr,"SimpleFeatureAdaptor_store not implemented\n");
  exit(1);

/*
  my ($self,@sf) = @_;

  if( scalar(@sf) == 0 ) {
    throw("Must call store with list of SimpleFeatures");
  }

  my $sth = $self->prepare
    ("INSERT INTO simple_feature (seq_region_id, seq_region_start, " .
                                 "seq_region_end, seq_region_strand, " .
                                 "display_label, analysis_id, score) " .
     "VALUES (?,?,?,?,?,?,?)");

  my $db = $self->db();
  my $analysis_adaptor = $db->get_AnalysisAdaptor();

 FEATURE: foreach my $sf ( @sf ) {

    if( !ref $sf || !$sf->isa("Bio::EnsEMBL::SimpleFeature") ) {
      throw("SimpleFeature must be an Ensembl SimpleFeature, " .
            "not a [".ref($sf)."]");
    }

    if($sf->is_stored($db)) {
      warning("SimpleFeature [".$sf->dbID."] is already stored" .
              " in this database.");
      next FEATURE;
    }

    if(!defined($sf->analysis)) {
      throw("An analysis must be attached to the features to be stored.");
    }

    #store the analysis if it has not been stored yet
    if(!$sf->analysis->is_stored($db)) {
      $analysis_adaptor->store($sf->analysis());
    }

    my $original = $sf;
    my $seq_region_id;
    ($sf, $seq_region_id) = $self->_pre_store($sf);

    $sth->bind_param(1,$seq_region_id,SQL_INTEGER);
    $sth->bind_param(2,$sf->start,SQL_INTEGER);
    $sth->bind_param(3,$sf->end,SQL_INTEGER);
    $sth->bind_param(4,$sf->strand,SQL_TINYINT);
    $sth->bind_param(5,$sf->display_label,SQL_VARCHAR);
    $sth->bind_param(6,$sf->analysis->dbID,SQL_INTEGER);
    $sth->bind_param(7,$sf->score,SQL_DOUBLE);

    $sth->execute();

    $original->dbID($sth->{'mysql_insertid'});
    $original->adaptor($self);
  }
*/
}


/*
=head2 _objs_from_sth

  Arg [1]    : hash reference $hashref
  Example    : none
  Description: PROTECTED implementation of superclass abstract method.
               creates SimpleFeatures from an executed DBI statement handle.
  Returntype : list reference to Bio::EnsEMBL::SimpleFeature objects
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut
*/
Vector *SimpleFeatureAdaptor_objectsFromStatementHandle(BaseFeatureAdaptor *bfa,
                                                        StatementHandle *sth,
                                                        AssemblyMapper *assMapper,
                                                        Slice *destSlice) {
  SliceAdaptor *sa     = DBAdaptor_getSliceAdaptor(bfa->dba);
  AnalysisAdaptor *aa  = DBAdaptor_getAnalysisAdaptor(bfa->dba);

  Vector *features = Vector_new();
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
    IDType simpleFeatureId  = row->getLongLongAt(row,0);
    IDType seqRegionId      = row->getLongLongAt(row,1);
    long seqRegionStart     = row->getLongAt(row,2);
    long seqRegionEnd       = row->getLongAt(row,3);
    int seqRegionStrand     = row->getIntAt(row,4);
    char *displayLabel      = row->getStringAt(row,5);
    IDType analysisId       = row->getLongLongAt(row,6);
    double score            = row->getDoubleAt(row,7);

    // get the analysis object
    Analysis *analysis = AnalysisAdaptor_fetchByDbID(aa, analysisId);

    if (! IDHash_contains(sliceHash, seqRegionId)) {
      IDHash_add(sliceHash, seqRegionId, SliceAdaptor_fetchBySeqRegionId(sa, seqRegionId, POS_UNDEF, POS_UNDEF, STRAND_UNDEF));
    }
    Slice *slice = IDHash_getValue(sliceHash, seqRegionId);

    Slice *sfSlice = slice;

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
      sfSlice = IDHash_getValue(sliceHash, seqRegionId);
    }

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
      sfSlice = destSlice;
    }

    SimpleFeature *sf = SimpleFeature_new();

    SimpleFeature_setSlice(sf, sfSlice);
    SimpleFeature_setStart(sf, seqRegionStart);
    SimpleFeature_setEnd(sf, seqRegionEnd);
    SimpleFeature_setStrand(sf, seqRegionStrand);
    SimpleFeature_setScore(sf, score);
    SimpleFeature_setAnalysis(sf,analysis);
    SimpleFeature_setAdaptor(sf, (BaseAdaptor *)bfa);
    SimpleFeature_setDbID(sf, simpleFeatureId);
    SimpleFeature_setDisplayLabel(sf, displayLabel);

    Vector_addElement(features, sf);
  }

  IDHash_free(sliceHash, NULL);
  return features;
}

/*
=head2 list_dbIDs

  Arg [1]    : none
  Example    : @feature_ids = @{$simple_feature_adaptor->list_dbIDs()};
  Description: Gets an array of internal ids for all simple features in the current db
  Arg[1]     : <optional> int. not 0 for the ids to be sorted by the seq_region.
  Returntype : list of ints
  Exceptions : none
  Caller     : ?
  Status     : Stable

=cut
*/
Vector *SimpleFeatureAdaptor_listDbIDs(SimpleFeatureAdaptor *sfa, int ordered) {

// NIY: Shouldn't really do this direct BaseAdaptor call but ...
  return BaseAdaptor_listDbIDs((BaseAdaptor *)sfa, "simple_feature", NULL, ordered);
}


