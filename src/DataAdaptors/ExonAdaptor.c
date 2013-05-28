#include "ExonAdaptor.h"
#include "AnalysisAdaptor.h"
#include "AssemblyMapperAdaptor.h"
#include "DBAdaptor.h"
#include "DBEntryAdaptor.h"
#include "StrUtil.h"
#include "BaseAdaptor.h"
#include "MysqlUtil.h"
#include "Exon.h"
#include "IDHash.h"
#include "Transcript.h"
#include "Slice.h"
#include "AssemblyMapper.h"
#include "ChainedAssemblyMapper.h"
#include "CoordSystemAdaptor.h"
#include "MetaCoordContainer.h"

#include "ExonAdaptor.h"
#include "TranscriptAdaptor.h"
#include "SliceAdaptor.h"

#include "StatementHandle.h"
#include "ResultRow.h"

#include "Error.h"
/*
=head1 DESCRIPTION

The ExonAdaptor is responsible for retrieving and storing Exon objects
from an Ensembl database.  Most of the ExonAdaptor functionality is
inherited from the B<Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor> class.
*/

static int ETMode = 0; // Hack to switch tables and final clause in fetchAllByTranscript
NameTableType ExonAdaptor_tableNamesStandard = {{"exon","e"},
                                                {NULL,NULL}};

NameTableType ExonAdaptor_tableNamesWithET   = {{"exon","e"},
                                                {"exon_transcript", "et"},
                                                {NULL,NULL}};


ExonAdaptor *ExonAdaptor_new(DBAdaptor *dba) {
  ExonAdaptor *ea;

  if ((ea = (ExonAdaptor *)calloc(1,sizeof(ExonAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for ExonAdaptor\n");
    return NULL;
  }
  BaseFeatureAdaptor_init((BaseFeatureAdaptor *)ea, dba, EXON_ADAPTOR);

  ea->getTables                  = ExonAdaptor_getTables;
  ea->getColumns                 = ExonAdaptor_getColumns;
  ea->store                      = ExonAdaptor_store;
  ea->objectsFromStatementHandle = ExonAdaptor_objectsFromStatementHandle;
  ea->finalClause                = ExonAdaptor_finalClause;

  return ea;
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
*/
NameTableType *ExonAdaptor_getTables() {
  if (ETMode == 0) {
    return &ExonAdaptor_tableNamesStandard;
  } else {
    return &ExonAdaptor_tableNamesWithET;
  }
}

/*
# _columns
#
#  Arg [1]    : none
#  Example    : none
#  Description: PROTECTED implementation of superclass abstract method
#               returns a list of columns to use for queries
#  Returntype : list of strings
#  Exceptions : none
#  Caller     : internal
*/
char *Exon_cols[] = {
                     "e.exon_id",
                     "e.seq_region_id",
                     "e.seq_region_start",
                     "e.seq_region_end",
                     "e.seq_region_strand",
                     "e.phase",
                     "e.end_phase",
                     "e.is_current",
                     "e.is_constitutive",
                     "e.stable_id",
                     "e.version",
                     "UNIX_TIMESTAMP(e.created_date)",
                     "UNIX_TIMESTAMP(e.modified_date)",
                     NULL };

char **ExonAdaptor_getColumns() {
  return Exon_cols;
}


/*
# _final_clause
#
#  Arg [1]    : none
#  Example    : none
#  Description: PROTECTED implementation of superclass abstract method
#               returns a default end for the SQL-query (ORDER BY)
#  Returntype : string
#  Exceptions : none
#  Caller     : internal
*/
char *ExonAdaptor_finalClause() {
  if (ETMode == 0) {
    return "";
  } else {
    return "ORDER BY et.transcript_id, et.rank";
  }
}

Vector *ExonAdaptor_fetchAll(ExonAdaptor *ea) {
  char constraint[1024];

// NIY: Maybe just a constant string
// Note: Doesn't ignore LRGs unlike genes and transcripts, so these bits of **** will get loaded
  sprintf(constraint, "e.is_current = 1");

  return ExonAdaptor_genericFetch(ea, constraint, NULL, NULL);
}

/*
=head2 fetch_by_stable_id

  Arg [1]    : string $stable_id
               the stable id of the exon to retrieve
  Example    : $exon = $exon_adaptor->fetch_by_stable_id('ENSE0000988221');
  Description: Retrieves an Exon from the database via its stable id
  Returntype : Bio::EnsEMBL::Exon in native coordinates.
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/
Exon *ExonAdaptor_fetchByStableId(ExonAdaptor *ea, char *stableId) {
  char constraint[1024];

  sprintf(constraint, "e.stable_id = '%s' AND e.is_current = 1", stableId);

  Vector *exons = ExonAdaptor_genericFetch(ea, constraint, NULL, NULL);
  Exon *exon = Vector_getElementAt(exons, 0);
// Need to free data for exons other than first one (same for gene and transcript in those adaptors)
  Vector_free(exons);

  return exon;
}


/*
=head2 fetch_all_versions_by_stable_id 

  Arg [1]     : String $stable_id 
                The stable ID of the exon to retrieve
  Example     : my $exon = $exon_adaptor->fetch_all_version_by_stable_id
                  ('ENSE00000309301');
  Description : Similar to fetch_by_stable_id, but retrieves all versions of an
                exon stored in the database.
  Returntype  : listref of Bio::EnsEMBL::Exon objects
  Exceptions  : if we cant get the gene in given coord system
  Caller      : general
  Status      : At Risk

=cut
*/
Vector *ExonAdaptor_fetchAllVersionsByStableId(ExonAdaptor *ea, char *stableId) {
  char constraint[1024];

  sprintf(constraint, "e.stable_id = '%s'", stableId);

  return ExonAdaptor_genericFetch(ea, constraint, NULL, NULL);
}


/*
=head2 fetch_all_by_Transcript

  Arg [1]    : Bio::EnsEMBL::Transcript $transcript
  Example    : none
  Description: Retrieves all Exons for the Transcript in 5-3 order
  Returntype : listref Bio::EnsEMBL::Exon on Transcript slice 
  Exceptions : throws if transcript has no slice
  Caller     : Transcript->get_all_Exons()
  Status     : Stable

=cut
*/
Vector *ExonAdaptor_fetchAllByTranscript(ExonAdaptor *ea, Transcript *transcript) {

  Slice *tSlice = Transcript_getSlice(transcript);
  Slice *slice;

  if (tSlice == NULL) {
    fprintf(stderr, "Transcript must have attached slice to retrieve exons.\n");
    exit(1);
  }

  // use a small slice the same size as the transcript
// No circular slice stuff
//  if ( !$tslice->is_circular() ) {
    SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(ea->dba);
    slice = SliceAdaptor_fetchByFeature(sa, (SeqFeature *)transcript, 0, 0);
//  } else {
//    # Circular.
//    $slice = $tslice;
//  }

  // Override the tables definition to provide an additional join to the
  // exon_transcript table.  For efficiency we cannot afford to have this
  // in as a left join every time.
  ETMode = 1; // Hacky flag to change behaviour of getTables and finalClause to add in the exon_transcript table

  char constraint[1024];
  sprintf(constraint, "et.transcript_id = "IDFMTSTR" AND e.exon_id = et.exon_id", Transcript_getDbID(transcript));

  // fetch all of the exons
  Vector *exons = ExonAdaptor_fetchAllBySliceConstraint(ea, slice, constraint, NULL);

  // un-override the table definition
  ETMode = 0;

  // remap exon coordinates if necessary
  if (EcoString_strcmp(Slice_getName(slice), Slice_getName(tSlice))) {
    Vector *out = Vector_new();
    int i;
    for (i=0; i<Vector_getNumElement(exons); i++) {
      Exon *ex = Vector_getElementAt(exons, i);
      Vector_addElement(out, Exon_transfer(ex, tSlice));
    }
// NIY Do I need to set a free func here - probably???
    Vector_free(exons);

    exons = out;
  }

  return exons;
}


/* NIY
=head2 store

  Arg [1]    : Bio::EnsEMBL::Exon $exon
               the exon to store in this database
  Example    : $exon_adaptor->store($exon);
  Description: Stores an exon in the database
  Returntype : none
  Exceptions : thrown if exon (or component exons) do not have a contig_id
               or if $exon->start, $exon->end, $exon->strand, or $exon->phase 
               are not defined or if $exon is not a Bio::EnsEMBL::Exon
  Caller     : general
  Status     : Stable

=cut
*/


IDType ExonAdaptor_store(ExonAdaptor *ea, Exon *exon) {
  if (exon == NULL) {
    fprintf(stderr, "feature is NULL in Exon_store\n");
    exit(1);
  }

  Class_assertType(CLASS_EXON, exon->objectType);

  DBAdaptor *db = ea->dba;

  if (Exon_isStored(exon, db)) {
    fprintf(stderr, "Exon ["IDFMTSTR"] is already stored in this database.\n", Exon_getDbID(exon) );
    return Exon_getDbID(exon);
  }

/* This check is odd - 0 would be OK for start or end if its on a slice which doesn't start at 1 in the seq region
   so this method is relying on being called after a transfer has been done in GeneAdaptor or TranscriptAdaptor.
   Doesn't seem correct. I'll do a cut down check
  if( ! $exon->start || ! $exon->end ||
      ! $exon->strand || ! defined $exon->phase ) {
    throw("Exon does not have all attributes to store");
  }
*/
  if (!Exon_getStrand(exon) || Exon_getPhase(exon) < -1 || Exon_getPhase(exon) > 2) {
    fprintf(stderr,"Exon does not have all attributes to store\n");
    exit(1);
  }

  // Default to is_current = 1 if this attribute is not set
  int isCurrent = Exon_getIsCurrent(exon);
/* Note in C isCurrent is initialised to 1 in Exon_new, so this check should be unnecessary 
  if ( !defined($is_current) ) { $is_current = 1 }
*/

  // Default to is_constitutive = 0 if this attribute is not set
  int isConstitutive = Exon_getIsConstitutive(exon);
/* Note in C isConstitutive will be 0 by default because Exon_new uses calloc to allocate the Exon
  if ( !defined($is_constitutive) ) { $is_constitutive = 0 }
*/

/* Moved up so have seqRegionId before making query string
   Not doing transfer in preStore so no new exon to worry about
  my $original = $exon;
  my $seq_region_id;
  ($exon, $seq_region_id) = $self->_pre_store($exon);
*/
  IDType seqRegionId = BaseFeatureAdaptor_preStore((BaseFeatureAdaptor *)ea, exon);

  char qStr[2048];
  strcpy(qStr,"INSERT into exon ( seq_region_id, seq_region_start, "
		       "seq_region_end, seq_region_strand, phase, "
		       "end_phase, is_current, is_constitutive");

  int version;
  if (Exon_getStableId(exon) != NULL) {
/*
    my $created = $self->db->dbc->from_seconds_to_date($exon->created_date());
    my $modified = $self->db->dbc->from_seconds_to_date($exon->modified_date());
*/
    version = Exon_getVersion(exon) <= 0 ? Exon_getVersion(exon) : 1;
    strcat(qStr,", stable_id, version, created_date, modified_date");
     
  }
  sprintf(qStr, "%s) VALUES ("IDFMTSTR", %ld, %ld, %d, %d, %d, %d, %d",
          qStr,
          (IDType)seqRegionId,
          Exon_getSeqRegionStart(exon),
          Exon_getSeqRegionEnd(exon),
          Exon_getSeqRegionStrand(exon),
          Exon_getPhase(exon),
          Exon_getEndPhase(exon),
          isCurrent,
          isConstitutive);
          
  if (Exon_getStableId(exon) != NULL) {
    sprintf(qStr, "%s, '%s', %d, FROM_UNIXTIME(%ld), FROM_UNIXTIME(%ld))",
            qStr,
            Exon_getStableId(exon),
            version,
            Exon_getCreated(exon),
            Exon_getModified(exon));
  } else {
    strcat(qStr,")");
  }

  StatementHandle *sth = ea->prepare((BaseAdaptor *)ea,qStr,strlen(qStr));

  IDType exonId = 0;

  // store the exon
  sth->execute(sth);

  exonId = sth->getInsertId(sth);

  // Now the supporting evidence
  SupportingFeatureAdaptor *esfAdaptor = DBAdaptor_getSupportingFeatureAdaptor(db);
  SupportingFeatureAdaptor_store(esfAdaptor, exonId, Exon_getAllSupportingFeatures(exon));

  // HISTORIC NOTE: This comment (the component exon bit) must be the last remnant of StickyExon code in
  // the API - aaahhhh. Thank goodness the little b***ers are gone.

  //
  // Finally, update the dbID and adaptor of the exon (and any component exons)
  // to point to the new database
  //
  Exon_setAdaptor(exon, (BaseAdaptor *)ea);
  Exon_setDbID(exon, exonId);

  return exonId;
}


/* NIY
=head2 remove

  Arg [1]    : Bio::EnsEMBL::Exon $exon
               the exon to remove from the database
  Example    : $exon_adaptor->remove($exon);
  Description: Removes an exon from the database.  This method is generally
               called by the TranscriptAdaptor::store method. Database
               integrity will not be maintained if this method is simply
               called on its own without taking into account transcripts which
               may refer to the exon being removed.
  Returntype : none
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub remove {
  my $self = shift;
  my $exon = shift;

  if(!ref($exon) || !$exon->isa('Bio::EnsEMBL::Exon')) {
    throw('Bio::EnsEMBL::Exon argument expected.');
  }

  if(!$exon->is_stored($self->db())) {
    warning("Cannot remove exon " .$exon->dbID.
            "Is not stored in this database.");
    return;
  }

  # sanity check: make sure nobdody tries to slip past a prediction exon
  # which inherits from exon but actually uses different tables
  if($exon->isa('Bio::EnsEMBL::PredictionExon')) {
    throw("ExonAdaptor can only remove Exons not PredictionExons.");
  }

  # Remove the supporting features of this exon

  my $prot_adp = $self->db->get_ProteinAlignFeatureAdaptor;
  my $dna_adp = $self->db->get_DnaAlignFeatureAdaptor;

  my $sth = $self->prepare("SELECT feature_type, feature_id  " .
                           "FROM supporting_feature " .            
			   "WHERE exon_id = ?");
  $sth->bind_param(1, $exon->dbID, SQL_INTEGER);
  $sth->execute();

  # statements to check for shared align_features
  my $sth1 = $self->prepare("SELECT count(*) FROM supporting_feature " .
			    "WHERE feature_type = ? AND feature_id = ?");
  my $sth2 = $self->prepare("SELECT count(*) " .
                            "FROM transcript_supporting_feature " .
			    "WHERE feature_type = ? AND feature_id = ?");

  SUPPORTING_FEATURE:
  while(my ($type, $feature_id) = $sth->fetchrow()){
    
    # only remove align_feature if this is the last reference to it
    $sth1->bind_param(1, $type, SQL_VARCHAR);
    $sth1->bind_param(2, $feature_id, SQL_INTEGER);
    $sth1->execute;
    $sth2->bind_param(1, $type, SQL_VARCHAR);
    $sth2->bind_param(2, $feature_id, SQL_INTEGER);
    $sth2->execute;
    my ($count1) = $sth1->fetchrow;
    my ($count2) = $sth2->fetchrow;
    if ($count1 + $count2 > 1) {
      #warn "shared feature, not removing $type|$feature_id\n";
      next SUPPORTING_FEATURE;
    }
    
    #warn "removing $type|$feature_id\n";
  
    if($type eq 'protein_align_feature'){
      my $f = $prot_adp->fetch_by_dbID($feature_id);
      $prot_adp->remove($f);
    }
    elsif($type eq 'dna_align_feature'){
      my $f = $dna_adp->fetch_by_dbID($feature_id);
      $dna_adp->remove($f);
    }
    else {
      warning("Unknown supporting feature type $type. Not removing feature.");
    }
  }
  $sth->finish();
  $sth1->finish();
  $sth2->finish();

  # delete the association to supporting features

  $sth = $self->prepare("DELETE FROM supporting_feature WHERE exon_id = ?");
  $sth->bind_param(1, $exon->dbID, SQL_INTEGER);
  $sth->execute();
  $sth->finish();


  # delete the exon

  $sth = $self->prepare( "DELETE FROM exon WHERE exon_id = ?" );
  $sth->bind_param(1, $exon->dbID, SQL_INTEGER);
  $sth->execute();
  $sth->finish();

  $exon->dbID(undef);
  $exon->adaptor(undef);

  return;
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
Vector *ExonAdaptor_listDbIDs(ExonAdaptor *ea, int ordered) {

// NIY: Shouldn't really do this direct BaseAdaptor call but ...
  return BaseAdaptor_listDbIDs((BaseAdaptor *)ea, "exon", NULL, ordered);
}


/*
=head2 list_stable_ids

  Arg [1]    : none
  Example    : @stable_exon_ids = @{$exon_adaptor->list_stable_dbIDs()};
  Description: Gets an array of stable ids for all exons in the current db
  Returntype : list of ints
  Exceptions : none
  Caller     : ?
  Status     : Stable

=cut
*/
Vector *ExonAdaptor_listStableIDs(ExonAdaptor *ea) {

// NIY: Shouldn't really do this direct BaseAdaptor call but ...
  return BaseAdaptor_listDbIDs((BaseAdaptor *)ea, "exon", "stable_id", 0);
}


/*
#_objs_from_sth
#
#  Arg [1]    : StatementHandle $sth
#  Example    : none 
#  Description: PROTECTED implementation of abstract superclass method.
#               responsible for the creation of Exons
#  Returntype : listref of Bio::EnsEMBL::Exons in target coordinate system
#  Exceptions : none
#  Caller     : internal
*/
Vector *ExonAdaptor_objectsFromStatementHandle(ExonAdaptor *ea,
                                               StatementHandle *sth,
                                               AssemblyMapper *assMapper,
                                               Slice *destSlice) {

  //
  // This code is ugly because an attempt has been made to remove as many
  // function calls as possible for speed purposes.  Thus many caches and
  // a fair bit of gymnastics is used.
  //

  SliceAdaptor *sa     = DBAdaptor_getSliceAdaptor(ea->dba);

  Vector *exons = Vector_new();
  IDHash *sliceHash = IDHash_new(IDHASH_SMALL);
/* Don't bother with these three - analysis is cached in its adaptor, and I can't believe speed
  my %sr_name_hash;
  my %sr_cs_hash;
*/



/* Basically Unused! - used but in a non sensible way (see TranscriptAdaptor for further comment)
  my $asm_cs;
  my $cmp_cs;
  my $asm_cs_vers;
  my $asm_cs_name;
  my $cmp_cs_vers;
  my $cmp_cs_name;

  if ($mapper) {
    $asm_cs      = $mapper->assembled_CoordSystem();
    $cmp_cs      = $mapper->component_CoordSystem();
    $asm_cs_name = $asm_cs->name();
    $asm_cs_vers = $asm_cs->version();
    $cmp_cs_name = $cmp_cs->name();
    $cmp_cs_vers = $cmp_cs->version();
  }
*/

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
    asma            = DBAdaptor_getAssemblyMapperAdaptor(ea->dba);
  }

// Note FEATURE label is here
//FEATURE: while ($sth->fetch())
  ResultRow *row;
// Extra parentheses to please mac compiler
  while ((row = sth->fetchRow(sth))) {
    IDType exonId       = row->getLongLongAt(row, 0);
    IDType seqRegionId  = row->getLongLongAt(row, 1);
    long seqRegionStart = row->getLongAt(row, 2);
    long seqRegionEnd   = row->getLongAt(row, 3);
    int seqRegionStrand = row->getIntAt(row, 4);
    int phase           = row->getIntAt(row, 5);
    int endPhase        = row->getIntAt(row, 6);
    int isCurrent       = row->getIntAt(row, 7);
    int isConstitutive  = row->getIntAt(row, 8);
    char *stableId      = row->getStringAt(row, 9);
    int version         = row->getIntAt(row, 10);
    int createdDate     = row->getIntAt(row, 11);
    int modifiedDate    = row->getIntAt(row, 12);

    // Not doing internal seq id stuff for now
//    #need to get the internal_seq_region, if present
//    $seq_region_id = $self->get_seq_region_id_internal($seq_region_id);
//    my $slice       = $slice_hash{ "ID:" . $seq_region_id };
//    my $dest_mapper = $mapper;
//    if ( !$slice ) {
//      $slice = $sa->fetch_by_seq_region_id($seq_region_id);
//      $slice_hash{ "ID:" . $seq_region_id } = $slice;
//      $sr_name_hash{$seq_region_id}         = $slice->seq_region_name();
//      $sr_cs_hash{$seq_region_id}           = $slice->coord_system();
//    }
//    my $sr_name = $sr_name_hash{$seq_region_id};
//    my $sr_cs   = $sr_cs_hash{$seq_region_id};

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
/*
      $asm_cs      = $dest_mapper->assembled_CoordSystem();
      $cmp_cs      = $dest_mapper->component_CoordSystem();
      $asm_cs_name = $asm_cs->name();
      $asm_cs_vers = $asm_cs->version();
      $cmp_cs_name = $cmp_cs->name();
      $cmp_cs_vers = $cmp_cs->version();
*/
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
        mrs = ChainedAssemblyMapper_map(assMapper, srName, seqRegionStart, seqRegionEnd, seqRegionStrand, srCs, 1, destSlice);
      } else {
        mrs = AssemblyMapper_fastMap(assMapper, srName, seqRegionStart, seqRegionEnd, seqRegionStrand, srCs, NULL);
      }

      // skip features that map to gaps or coord system boundaries
      //next FEATURE if (!defined($seq_region_id));
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
    // If a destination slice was provided convert the coords.
    //  
    if (destSlice != NULL) {
      if (destSliceStrand == 1) {
        // Positive strand.
        seqRegionStart = seqRegionStart - destSliceStart + 1;
        seqRegionEnd   = seqRegionEnd - destSliceStart + 1;

        if (0) {
/*
	if ( ( $seq_region_end > $dest_slice_start || $seq_region_end < 0 || ( $dest_slice_start > $dest_slice_end
                 && $seq_region_end < 0 ) )  && $dest_slice->is_circular() ) {
          # Handle circular chromosomes.

          if ( $seq_region_start > $seq_region_end ) {
            # Looking at a feature overlapping the chromsome origin.

            if ( $seq_region_end > $dest_slice_start ) {
              # Looking at the region in the beginning of the
              # chromosome.
              $seq_region_start -= $dest_slice->seq_region_length();
            }

            if ( $seq_region_end < 0 ) {
              $seq_region_end += $dest_slice->seq_region_length();
            }

          } else {
            if (    $dest_slice_start > $dest_slice_end
                 && $seq_region_end < 0 )
            {
              # Looking at the region overlapping the chromosome
              # origin and a feature which is at the beginning of the
              # chromosome.
              $seq_region_start += $dest_slice->seq_region_length();
              $seq_region_end   += $dest_slice->seq_region_length();
            }
          }
*/
        }
      } else {
        // Negative strand.
        if (0) {
/*
        if ( $seq_region_start > $seq_region_end && $dest_slice->is_circular() )
          # Handle circular chromosomes.

          if ( $dest_slice_start > $dest_slice_end ) {
            my $tmp_seq_region_start = $seq_region_start;
            $seq_region_start = $dest_slice_end - $seq_region_end + 1;
            $seq_region_end =
              $dest_slice_end +
              $dest_slice->seq_region_length() -
              $tmp_seq_region_start + 1;
          } else {

            if ( $seq_region_end > $dest_slice_start ) {
              # Looking at the region in the beginning of the
              # chromosome.
              $seq_region_start = $dest_slice_end - $seq_region_end + 1;
              $seq_region_end =
                $seq_region_end -
                $dest_slice->seq_region_length() -
                $dest_slice_start + 1;
            } else {
              my $tmp_seq_region_start = $seq_region_start;
              $seq_region_start =
                $dest_slice_end -
                $seq_region_end -
                $dest_slice->seq_region_length() + 1;
              $seq_region_end =
                $dest_slice_end - $tmp_seq_region_start + 1;
            }

          }
*/

        } else {
          // Non-circular chromosome - sanity
          long tmpSeqRegionStart = seqRegionStart;
          seqRegionStart = destSliceEnd - seqRegionEnd + 1;
          seqRegionEnd   = destSliceEnd - tmpSeqRegionStart + 1;
        }

        seqRegionStrand = -seqRegionStrand;
      }

      // Throw away features off the end of the requested slice or on
      // different seq_region.
// Perl used 'ne' for comparison of ids but these are ints so that's not really very efficient
      if (seqRegionEnd < 1 ||
          seqRegionStart > destSliceLength ||
          (destSliceSrId != seqRegionId)) {
// Any freeing to do - don't think so??
        //next FEATURE;
        continue;
      }
      exonSlice = destSlice;
    }

    // Finally, create the new exon.
    Exon *exon = Exon_new();
    Exon_setStart          (exon, seqRegionStart);
    Exon_setEnd            (exon, seqRegionEnd);
    Exon_setStrand         (exon, seqRegionStrand);
    Exon_setAdaptor        (exon, (BaseAdaptor *)ea);
    Exon_setSlice          (exon, exonSlice);
    Exon_setDbID           (exon, exonId);
    Exon_setStableId       (exon, stableId);
    Exon_setVersion        (exon, version);
    Exon_setCreated        (exon, createdDate);
    Exon_setModified       (exon, modifiedDate);
    Exon_setPhase          (exon, phase);
    Exon_setEndPhase       (exon, endPhase);
    Exon_setIsCurrent      (exon, isCurrent);
    Exon_setIsConstitutive (exon, isCurrent);
 
    Vector_addElement(exons, exon);
  }

  // Don't free slices because they might be being used????
  IDHash_free(sliceHash, NULL);

  return exons;
}

/* Don't bother
=head1 DEPRECATED METHODS

=cut


=head2 get_stable_entry_info

  Description: DEPRECATED. This method is no longer necessary.  Exons are
               always fetched with their stable identifiers (if they exist) and
               no lazy loading is necessary.

=cut

sub get_stable_entry_info {
  my ($self,$exon) = @_;

  deprecated( "This method call shouldnt be necessary" );

  if( !$exon || !ref $exon || !$exon->isa('Bio::EnsEMBL::Exon') ) {
     $self->throw("Needs a exon object, not a $exon");
  }
  if(!$exon->dbID){
    #$self->throw("can't fetch stable info with no dbID");
    return;
  }

  my $created_date = $self->db->dbc->from_date_to_seconds("created_date");
  my $modified_date = $self->db->dbc->from_date_to_seconds("modified_date");
  my $sth = $self->prepare("SELECT stable_id, " . $created_date . ",
                                   " . $modified_date . ", version 
                            FROM   exon
                            WHERE  exon_id = ");

  $sth->bind_param(1, $exon->dbID, SQL_INTEGER);
  $sth->execute();

  # my @array = $sth->fetchrow_array();
  if( my $aref = $sth->fetchrow_arrayref() ) {
    $exon->{'_stable_id'} = $aref->[0];
    $exon->{'_created'}   = $aref->[1];
    $exon->{'_modified'}  = $aref->[2];
    $exon->{'_version'}   = $aref->[3];
  }

  return 1;
}


=head2 fetch_all_by_gene_id

  Description: DEPRECATED. This method should not be needed - Exons can
               be fetched by Transcript.

=cut

sub fetch_all_by_gene_id {
  my ( $self, $gene_id ) = @_;
  my %exons;
  my $hashRef;
  my ( $currentId, $currentTranscript );

  deprecated( "Hopefully this method is not needed any more. Exons should be fetched by Transcript" );

  if( !$gene_id ) {
      $self->throw("Gene dbID not defined");
  }
  
  $self->{rchash} = {};
  
  my $query = qq {
    SELECT 
      STRAIGHT_JOIN 
	e.exon_id
      , e.contig_id
      , e.contig_start
      , e.contig_end
      , e.contig_strand
      , e.phase
      , e.end_phase
      , e.sticky_rank
    FROM transcript t
      , exon_transcript et
      , exon e
    WHERE t.gene_id = ?
      AND et.transcript_id = t.transcript_id
      AND e.exon_id = et.exon_id
    ORDER BY t.transcript_id,e.exon_id
      , e.sticky_rank DESC
  };

  my $sth = $self->prepare( $query );
  $sth->bind_param(1,$gene_id,SQL_INTEGER);
  $sth->execute();

  while( $hashRef = $sth->fetchrow_hashref() ) {
    if( ! exists $exons{ $hashRef->{exon_id} } ) {

      my $exon = $self->_exon_from_sth( $sth, $hashRef );

      $exons{$exon->dbID} = $exon;
    }
  }
  delete $self->{rchash};
  
  my @out = ();

  push @out, values %exons;

  return \@out;
}

*/
