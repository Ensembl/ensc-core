/*
=head1 DESCRIPTION

This is an adaptor responsible for the retrieval and storage of
DnaDnaAlignFeatures from the database. This adaptor inherits most of its
functionality from the BaseAlignFeatureAdaptor and BaseFeatureAdaptor
superclasses.

=head1 METHODS
*/
#include "DNAAlignFeatureAdaptor.h"

#include "AnalysisAdaptor.h"
#include "DNAAlignFeature.h"
#include "SliceAdaptor.h"
#include "ChainedAssemblyMapper.h"

NameTableType DNAAlignFeatureAdaptor_tableNames = {{"dna_align_feature","daf"},
                                                   {"external_db", "exdb"},
                                                   {NULL, NULL}};


NameTableType DNAAlignFeatureAdaptor_leftJoins = {{"external_db","exdb.external_db_id = daf.external_db_id"},
                                                  {NULL,NULL}};

DNAAlignFeatureAdaptor *DNAAlignFeatureAdaptor_new(DBAdaptor *dba) {
  DNAAlignFeatureAdaptor *dafa;

  if ((dafa = (DNAAlignFeatureAdaptor *)calloc(1,sizeof(DNAAlignFeatureAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for DNAAlignFeatureAdaptor\n");
    return NULL;
  }
  BaseFeatureAdaptor_init((BaseFeatureAdaptor *)dafa, dba, DNAALIGNFEATURE_ADAPTOR);

  dafa->getTables                  = DNAAlignFeatureAdaptor_getTables;
  dafa->getColumns                 = DNAAlignFeatureAdaptor_getColumns;
  dafa->store                      = DNAAlignFeatureAdaptor_store;
  dafa->objectsFromStatementHandle = DNAAlignFeatureAdaptor_objectsFromStatementHandle;
  dafa->leftJoin                   = DNAAlignFeatureAdaptor_leftJoin;

  return dafa;
}


/*
=head2 _tables

  Args       : none
  Example    : @tabs = $self->_tables
  Description: PROTECTED implementation of the abstract method inherited from
               BaseFeatureAdaptor.  Returns list of [tablename, alias] pairs
  Returntype : list of listrefs of strings
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor::generic_fetch
  Status     : Stable

=cut
*/
NameTableType *DNAAlignFeatureAdaptor_getTables(void) {
  return &DNAAlignFeatureAdaptor_tableNames;
}


NameTableType *DNAAlignFeatureAdaptor_leftJoin(void) {
  return &DNAAlignFeatureAdaptor_leftJoins;
}


/*
=head2 _columns

  Args       : none
  Example    : @columns = $self->_columns
  Description: PROTECTED implementation of abstract superclass method.  
               Returns a list of columns that are needed for object creation.
  Returntype : list of strings
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor::generic_fetch
  Status     : Stable

=cut
*/
char *DNAAlign_cols[] = {
            "daf.dna_align_feature_id",
            "daf.seq_region_id",
            "daf.analysis_id",
            "daf.seq_region_start",
            "daf.seq_region_end",
            "daf.seq_region_strand",
            "daf.hit_start",
            "daf.hit_end",
            "daf.hit_name",
            "daf.hit_strand",
            "daf.cigar_line",
            "daf.evalue",
            "daf.perc_ident",
            "daf.score",
            "daf.external_db_id",
            "daf.hcoverage",
	    "daf.external_data",
	    "daf.pair_dna_align_feature_id",
	    "exdb.db_name",
	    "exdb.db_display_name",
            NULL };

char **DNAAlignFeatureAdaptor_getColumns(void) {
  return DNAAlign_cols;
}


/*
=head2 store

  Arg [1]    : list of Bio::EnsEMBL::DnaAlignFeatures @feats
               the features to store in the database
  Example    : $dna_align_feature_adaptor->store(@features);
  Description: Stores a list of DnaAlignFeatures in the database
  Returntype : none
  Exceptions : throw if any of the provided features cannot be stored
               which may occur if:
                 * The feature does not have an associate Slice
                 * The feature does not have an associated analysis
                 * The Slice the feature is associated with is on a seq_region
                   unknown to this database
               A warning is given if:
                 * The feature has already been stored in this db
  Caller     : Pipeline
  Status     : Stable

=cut
*/
int DNAAlignFeatureAdaptor_store(BaseFeatureAdaptor *bfa, Vector *features) {

  if (features == NULL || Vector_getNumElement(features) == 0) {
    fprintf(stderr,"Must call store with features\n");
    exit(1);
  }

  NameTableType *tables = bfa->getTables();
  char *tableName = (*tables)[0][NAME];

  DBAdaptor *db                    = bfa->dba;
  AnalysisAdaptor *analysisAdaptor = DBAdaptor_getAnalysisAdaptor(db);

// Not very happy with this way of doing this - the undef values cause this to be ugly. Maybe should use bind parameters for storing
  char qStr[1024];
  sprintf(qStr,
    "INSERT INTO %s (seq_region_id,"
                    "seq_region_start,"
                    "seq_region_end,"
                    "seq_region_strand,"
                    "hit_start,"
                    "hit_end,"
                    "hit_strand,"
                    "hit_name,"
                    "cigar_line,"
                    "analysis_id,"
                    "score,"
                    "evalue,"
                    "perc_ident,"
                    "external_db_id,"
                    "hcoverage,"
                    "pair_dna_align_feature_id) "
                    "VALUES (%" IDFMTSTR ",%%d,%%d,%%d,%%d,%%d,%%d,'%%s','%%s',%"IDFMTSTR
                             ",%%s,%%s,%%s,%%s,%%s,%%s)", tableName);
                            // changed to strings to allow for nulls ",%%f,%%f,%%f,%" IDFMTSTR ",%%f,%" IDFMTSTR ")", tableName);

  StatementHandle *sth = bfa->prepare((BaseAdaptor *)bfa, qStr,strlen(qStr));

  int i;
  for (i=0; i<Vector_getNumElement(features); i++) {
    char fixedCigar[1024];

    char qStr[2048];

    DNAAlignFeature *feat = Vector_getElementAt(features, i);

    if (feat == NULL) {
      fprintf(stderr, "feature is NULL in DNAAlignFeature_store\n");
      exit(1);
    }

    Class_assertType(CLASS_DNADNAALIGNFEATURE, feat->objectType);

    if (DNAAlignFeature_isStored(feat, db)) {
      fprintf(stderr, "DNAAlignFeature ["IDFMTSTR"] is already stored in this database.\n", DNAAlignFeature_getDbID(feat) );
      continue;
    }

    BaseFeatureAdaptor_checkStartEndStrand(bfa,
                                           DNAAlignFeature_getHitStart(feat),
                                           DNAAlignFeature_getHitEnd(feat),
                                           DNAAlignFeature_getHitStrand(feat),
                                           NULL);

    char *cigarString = DNAAlignFeature_getCigarString(feat);

    if (cigarString == NULL) {
      sprintf(fixedCigar, "%ldM", DNAAlignFeature_getLength(feat));
      cigarString = fixedCigar;
      fprintf(stderr, "DNAAlignFeature does not define a cigar_string.\n"
                      "Assuming ungapped block with cigar_line = %s.\n", cigarString);
    }

    if (DNAAlignFeature_getHitSeqName(feat) == NULL) {
      fprintf(stderr, "DNAAlignFeature must define an hseqname.\n");
      exit(1);
    }

    Analysis *analysis = DNAAlignFeature_getAnalysis(feat);
    if (analysis == NULL) {
      fprintf(stderr,"An analysis must be attached to the features to be stored.\n");
      exit(1);
    }

    // store the analysis if it has not been stored yet
    if (Analysis_isStored(analysis, db)) {
      AnalysisAdaptor_store(analysisAdaptor, analysis);
    }

// Note no feature transfer in preStore currently so no new feature
// I can't see an obvious absolute requirement for the transfer to
// happen, so I'd really really really prefer NOT to do it, as it
// will cause all sorts of pain freeing up the temporary stuff
    //my $original = $feat;
    // ( $feat, $seq_region_id ) = $self->_pre_store($feat);
    IDType seqRegionId = BaseFeatureAdaptor_preStore(bfa, feat);

    char scoreString[1024];
    DNAAlignFeature_getScore(feat) != FLOAT_UNDEF ? sprintf(scoreString,"%f",DNAAlignFeature_getScore(feat)) : sprintf(scoreString,"NULL");

    char pValueString[1024];
    DNAAlignFeature_getpValue(feat) != FLOAT_UNDEF ? sprintf(pValueString,"%f",DNAAlignFeature_getpValue(feat)) : sprintf(pValueString,"NULL");

    char percIdString[1024];
    DNAAlignFeature_getPercId(feat) != FLOAT_UNDEF ? sprintf(percIdString,"%f",DNAAlignFeature_getPercId(feat)) : sprintf(percIdString,"NULL");

    char hCoverageString[1024];
    DNAAlignFeature_gethCoverage(feat) != FLOAT_UNDEF ? sprintf(hCoverageString,"%f",DNAAlignFeature_gethCoverage(feat)) : sprintf(hCoverageString,"NULL");

    char externalDbIDString[1024];
    DNAAlignFeature_getExternalDbID(feat) != 0 ? sprintf(externalDbIDString,IDFMTSTR,DNAAlignFeature_getExternalDbID(feat)) : sprintf(externalDbIDString,"NULL");

    char pairDNAAlignIDString[1024];
    DNAAlignFeature_getPairDNAAlignFeatureId(feat) != 0 ? sprintf(pairDNAAlignIDString,IDFMTSTR,DNAAlignFeature_getPairDNAAlignFeatureId(feat)) : sprintf(pairDNAAlignIDString,"NULL");

// Note using SeqRegionStart etc here rather than Start - should have same effect as perl's transfer
    
    sth->execute(sth, (IDType)seqRegionId,
                      DNAAlignFeature_getSeqRegionStart(feat),
                      DNAAlignFeature_getSeqRegionEnd(feat),
                      DNAAlignFeature_getSeqRegionStrand(feat),
                      DNAAlignFeature_getHitStart(feat),
                      DNAAlignFeature_getHitEnd(feat),
                      DNAAlignFeature_getHitStrand(feat),
                      DNAAlignFeature_getHitSeqName(feat),
                      cigarString,
                      (IDType)Analysis_getDbID(analysis),
                      scoreString,
                      pValueString,
                      percIdString,
                      externalDbIDString,
                      hCoverageString,
                      pairDNAAlignIDString);
       
/*
                      DNAAlignFeature_getScore(feat),
                      DNAAlignFeature_getpValue(feat),
                      DNAAlignFeature_getPercId(feat),
                      (IDType)DNAAlignFeature_getExternalDbID(feat),
                      DNAAlignFeature_gethCoverage(feat),
                      (IDType)DNAAlignFeature_getPairDNAAlignFeatureId(feat));
*/


    DNAAlignFeature_setDbID(feat,sth->getInsertId(sth));
    DNAAlignFeature_setAdaptor(feat, (BaseAdaptor *)bfa);
  }

  sth->finish(sth);
}


/* Huhhh?? what is this
sub save {
  my ($self, $features) = @_;

  my @feats = @$features;
  throw("Must call store with features") if( scalar(@feats) == 0 );

  my @tabs = $self->_tables;
  my ($tablename) = @{$tabs[0]};

  my $db = $self->db();
  my $analysis_adaptor = $db->get_AnalysisAdaptor();

  my $sql = qq{INSERT INTO $tablename (seq_region_id, seq_region_start, seq_region_end, seq_region_strand, hit_start, hit_end, hit_strand, hit_name, cigar_line, analysis_id, score, evalue, perc_ident, external_db_id, hcoverage, pair_dna_align_feature_id, external_data) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)};

  my %analyses = ();

  my $sth = $self->prepare($sql);
     
 FEATURE: foreach my $feat ( @feats ) {
    if( !ref $feat || !$feat->isa("Bio::EnsEMBL::DnaDnaAlignFeature") ) {
      throw("feature must be a Bio::EnsEMBL::DnaDnaAlignFeature,"
            . " not a [".ref($feat)."].");
    }

    if($feat->is_stored($db)) {
      warning("DnaDnaAlignFeature [".$feat->dbID."] is already stored" .
              " in this database.");
      next FEATURE;
    }

    my $hstart  = $feat->hstart || 0; # defined $feat->hstart  ? $feat->hstart : $feat->start ;
    my $hend    = $feat->hend   || 0; # defined $feat->hend    ? $feat->hend : $feat->end;
    my $hstrand = $feat->hstrand|| 0; # defined $feat->hstrand ? $feat->hstrand : $feat->strand;
    if( $hstart && $hend ) {
      if($hend < $hstart) {
        throw("Invalid Feature start/end [$hstart/$hend]. Start must be less than or equal to end.");
      }
    }
    my $cigar_string = $feat->cigar_string();
    if(!$cigar_string) {
      $cigar_string = $feat->length() . 'M';
      warning("DnaDnaAlignFeature does not define a cigar_string.\n" .
              "Assuming ungapped block with cigar_line=$cigar_string .");
    }

    my $hseqname = $feat->hseqname();
    if(!$hseqname) {
      throw("DnaDnaAlignFeature must define an hseqname.");
    }

    if(!defined($feat->analysis)) {
      throw("An analysis must be attached to the features to be stored.");
    }

    #store the analysis if it has not been stored yet
    if(!$feat->analysis->is_stored($db)) {
      $analysis_adaptor->store($feat->analysis());
    }

    $analyses{ $feat->analysis->dbID }++;

    my $original = $feat;
    my $seq_region_id;
    ($feat, $seq_region_id) = $self->_pre_store_userdata($feat);

    my $extra_data = $feat->extra_data ? $self->dump_data($feat->extra_data) : '';

    $sth->bind_param(1,$seq_region_id,SQL_INTEGER);
    $sth->bind_param(2,$feat->start,SQL_INTEGER);
    $sth->bind_param(3,$feat->end,SQL_INTEGER);
    $sth->bind_param(4,$feat->strand,SQL_TINYINT);
    $sth->bind_param(5,$hstart,SQL_INTEGER);
    $sth->bind_param(6,$hend,SQL_INTEGER);
    $sth->bind_param(7,$hstrand,SQL_TINYINT);
    $sth->bind_param(8,$hseqname,SQL_VARCHAR);
    $sth->bind_param(9,$cigar_string,SQL_LONGVARCHAR);
    $sth->bind_param(10,$feat->analysis->dbID,SQL_INTEGER);
    $sth->bind_param(11,$feat->score,SQL_DOUBLE);
#    $sth->bind_param(11,$feat->score); # if the above statement does not work it means you need to upgrade DBD::mysql, meantime you can replace it with this line
    $sth->bind_param(12,$feat->p_value,SQL_DOUBLE);
    $sth->bind_param(13,$feat->percent_id,SQL_FLOAT);
    $sth->bind_param(14,$feat->external_db_id,SQL_INTEGER);
    $sth->bind_param(15,$feat->hcoverage,SQL_DOUBLE);
    $sth->bind_param(16,$feat->pair_dna_align_feature_id,SQL_INTEGER);
    $sth->bind_param(17,$extra_data,SQL_LONGVARCHAR);


    $sth->execute();
    $original->dbID($sth->{'mysql_insertid'});
    $original->adaptor($self);
  }

  $sth->finish();

## js5 hack to update meta_coord table... 
  if( keys %analyses ) {

    my $sth = $self->prepare( 'select sr.coord_system_id, max(daf.seq_region_end-daf.seq_region_start) from seq_region as sr, dna_align_feature as daf where daf.seq_region_id=sr.seq_region_id and analysis_id in ('.join(',',keys %analyses).') group by coord_system_id' );
    $sth->execute;

    foreach( @{ $sth->fetchall_arrayref } ) {
      my $sth2 = $self->prepare( qq(insert ignore into meta_coord values("dna_align_feature",$_->[0],$_->[1])) );
      $sth2->execute;
      $sth2->finish;

      $sth2 = $self->prepare( qq(update meta_coord set max_length = $_->[1] where coord_system_id = $_->[0] and table_name="dna_align_feature" and max_length < $_->[1]) );
      $sth2->execute;
      $sth2->finish;
    }

    $sth->finish;
  }

}
*/


/*
=head2 _objs_from_sth

  Arg [1]    : DBI statement handle $sth
               an exectuted DBI statement handle generated by selecting 
               the columns specified by _columns() from the table specified 
               by _table()
  Example    : @dna_dna_align_feats = $self->_obj_from_hashref
  Description: PROTECTED implementation of superclass abstract method. 
               Creates DnaDnaAlignFeature objects from a DBI hashref
  Returntype : listref of Bio::EnsEMBL::DnaDnaAlignFeatures
  Exceptions : none
  Caller     : Bio::EnsEMBL::BaseFeatureAdaptor::generic_fetch
  Status     : Stable

=cut
*/
Vector *DNAAlignFeatureAdaptor_objectsFromStatementHandle(BaseFeatureAdaptor *bfa,
                                                       StatementHandle *sth,
                                                       AssemblyMapper *assMapper,
                                                       Slice *destSlice) {
/* Not doing this specialness
  # In case of userdata we need the features on the dest_slice.  In case
  # of get_all_supporting_features dest_slice is not provided.
  my $sa = (   $dest_slice
             ? $dest_slice->adaptor()
             : $self->db()->get_SliceAdaptor() );
*/
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
  while ((row = sth->fetchRow(sth))) {
    IDType dnaAlignFeatureId     = row->getLongLongAt(row,0);
    IDType seqRegionId           = row->getLongLongAt(row,1);
    IDType analysisId            = row->getLongLongAt(row,2);
    long seqRegionStart          = row->getLongAt(row,3);
    long seqRegionEnd            = row->getLongAt(row,4);
    int seqRegionStrand          = row->getIntAt(row,5);
    int hitStart                 = row->getIntAt(row,6);
    int hitEnd                   = row->getIntAt(row,7);
    char *hitName                = row->getStringAt(row,8);
    int hitStrand                = row->getIntAt(row,9);
    char *cigarLine              = row->getStringAt(row,10);
    double eValue                = row->getDoubleAt(row,11);
    double percIdent             = row->getDoubleAt(row,12);
    double score                 = row->getDoubleAt(row,13);
    IDType externalDbId          = row->getLongLongAt(row,14);
    double hCoverage             = row->getDoubleAt(row,15);
    char *extraData              = row->getStringAt(row,16);
    IDType pairDnaAlignFeatureId = row->getLongLongAt(row,17);
    char *externalDbName         = row->getStringAt(row,18);
    char *externalDbDisplayName  = row->getStringAt(row,19);

    // get the analysis object
    Analysis *analysis = AnalysisAdaptor_fetchByDbID(aa, analysisId);

    if (! IDHash_contains(sliceHash, seqRegionId)) {
      IDHash_add(sliceHash, seqRegionId, SliceAdaptor_fetchBySeqRegionId(sa, seqRegionId, POS_UNDEF, POS_UNDEF, STRAND_UNDEF));
    }
    Slice *slice = IDHash_getValue(sliceHash, seqRegionId);

    Slice *dafSlice = slice;

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
        mrs = ChainedAssemblyMapper_map(assMapper, srName, seqRegionStart, seqRegionEnd, seqRegionStrand, srCs, 1, destSlice);
      } else {
        mrs = AssemblyMapper_fastMap(assMapper, srName, seqRegionStart, seqRegionEnd, seqRegionStrand, srCs, NULL);
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
      dafSlice = IDHash_getValue(sliceHash, seqRegionId);
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
      dafSlice = destSlice;
    }


/* Wierd evil eval stuff - not doing
    # Inlining the following in the hash causes major issues with 5.16 and messes up the hash 
    my $evalled_extra_data = $extra_data ? $self->get_dumped_data($extra_data) : '';
*/

    // Finally, create the new DnaAlignFeature.
    DNAAlignFeature *daf = DNAAlignFeature_new();

    DNAAlignFeature_setSlice(daf,dafSlice);
    DNAAlignFeature_setStart(daf, seqRegionStart);
    DNAAlignFeature_setEnd(daf, seqRegionEnd);
    DNAAlignFeature_setStrand(daf, seqRegionStrand);

    DNAAlignFeature_setHitSeqName(daf,hitName);
    DNAAlignFeature_setHitStart(daf, hitStart);
    DNAAlignFeature_setHitEnd(daf, hitEnd);
    DNAAlignFeature_setHitStrand(daf, hitStrand);

    DNAAlignFeature_setScore(daf, score);
    DNAAlignFeature_setpValue(daf, eValue);
    DNAAlignFeature_setPercId(daf, percIdent);

    DNAAlignFeature_setCigarString(daf, cigarLine);

    DNAAlignFeature_setAnalysis(daf,analysis);
    DNAAlignFeature_setAdaptor(daf, (BaseAdaptor *)bfa);
    DNAAlignFeature_setDbID(daf, dnaAlignFeatureId);

    DNAAlignFeature_setExternalDbID(daf, externalDbId);
    DNAAlignFeature_sethCoverage(daf, hCoverage);
    // Unevaled in C           'extra_data'     => $evalled_extra_data,
    DNAAlignFeature_setExtraData(daf, extraData);
    DNAAlignFeature_setDbName(daf, externalDbName);
    DNAAlignFeature_setDbDisplayName(daf, externalDbDisplayName);

    DNAAlignFeature_setPairDNAAlignFeatureId(daf, pairDnaAlignFeatureId);

    Vector_addElement(features, daf); 
  } 

  IDHash_free(sliceHash, NULL);
  return features;
}

/*
=head2 list_dbIDs

  Arg [1]    : none
  Example    : @feature_ids = @{$dna_align_feature_adaptor->list_dbIDs()};
  Description: Gets an array of internal ids for all dna align features in 
               the current db
  Arg[1]     : <optional> int. not 0 for the ids to be sorted by the seq_region.
  Returntype : list of ints
  Exceptions : none
  Caller     : ?
  Status     : Stable

=cut
*/
Vector *DNAAlignFeatureAdaptor_listDbIDs(DNAAlignFeatureAdaptor *dafa, int ordered) {

// NIY: Shouldn't really do this direct BaseAdaptor call but ...
  return BaseAdaptor_listDbIDs((BaseAdaptor *)dafa, "dna_align_feature", NULL, ordered);
}

