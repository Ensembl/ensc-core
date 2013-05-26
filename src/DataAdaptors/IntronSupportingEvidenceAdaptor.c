/*
=head1 NAME

Bio::EnsEMBL::DBSQL::IntronSupportingEvidenceAdaptor

=head1 SYNOPSIS

  my $isea = $dba->get_IntronSupportingEvidenceAdaptor();
  my $ise = $isea->fetch_by_dbID(1);
  my $ise_array = $dfa->fetch_all();

*/

#include "IntronSupportingEvidenceAdaptor.h"

#include "DNAAlignFeatureAdaptor.h"
#include "ProteinAlignFeatureAdaptor.h"
#include "BaseAlignFeature.h"
#include "DBAdaptor.h"
#include "BaseFeatureAdaptor.h"

#include "AnalysisAdaptor.h"
#include "SliceAdaptor.h"

#include "AssemblyMapper.h"
#include "ChainedAssemblyMapper.h"

#include <string.h>

NameTableType IntronSupportingEvidenceAdaptor_tableNames = {{"intron_supporting_evidence","ise"},
                                                            {NULL,NULL}};

IntronSupportingEvidenceAdaptor *IntronSupportingEvidenceAdaptor_new(DBAdaptor *dba) {
  IntronSupportingEvidenceAdaptor *isea;

  if ((isea = (IntronSupportingEvidenceAdaptor *)calloc(1,sizeof(IntronSupportingEvidenceAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for IntronSupportingEvidenceAdaptor\n");
    return NULL;
  }
  BaseFeatureAdaptor_init((BaseFeatureAdaptor *)isea, dba, INTRONSUPPORTINGEVIDENCE_ADAPTOR);

  isea->getTables                  = IntronSupportingEvidenceAdaptor_getTables;
  isea->getColumns                 = IntronSupportingEvidenceAdaptor_getColumns;
  isea->store                      = IntronSupportingEvidenceAdaptor_store;
  isea->objectsFromStatementHandle = IntronSupportingEvidenceAdaptor_objectsFromStatementHandle;

  return isea;
}


NameTableType *IntronSupportingEvidenceAdaptor_getTables() {
  return &IntronSupportingEvidenceAdaptor_tableNames;
}

char *IntronSupportingEvidence_cols[] = { 
                                         "ise.intron_supporting_evidence_id",
                                         "ise.analysis_id",
                                         "ise.seq_region_id",
                                         "ise.seq_region_start",
                                         "ise.seq_region_end",
                                         "ise.seq_region_strand",
                                         "ise.hit_name",
                                         "ise.score",
                                         "ise.score_type",
                                         "ise.is_splice_canonical",
                                         NULL };

char **IntronSupportingEvidenceAdaptor_getColumns() {
  return IntronSupportingEvidence_cols;
}

/*
=head2 list_linked_transcript_ids

  Arg[1]      : Bio::EnsEMBL::IntronSupportingEvidence Evidence to search with
  Example     : my $transcript_ids = @{$isea->list_linked_transcript_ids($ise)};
  Description : Uses the given IntronSupportingEvidence to find all linked
                transcript ids 
  Returntype  : ArrayRef[Integer] of transcript_id
  Exceptions  : Thrown if arguments are not as stated and for DB errors 

=cut
*/

Vector *IntronSupportingEvidenceAdaptor_listLinkedTranscriptIds(IntronSupportingEvidenceAdaptor *isea, IntronSupportingEvidence *ise) {
  char qStr[1024];

  sprintf(qStr,"SELECT transcript_id from transcript_intron_supporting_evidence "
               "WHERE intron_supporting_evidence_id = "IDFMTSTR, IntronSupportingEvidence_getDbID(ise));

  StatementHandle *sth = isea->prepare((BaseAdaptor *)isea,qStr,strlen(qStr));
  sth->execute(sth);

  Vector *idVec = Vector_new();
  ResultRow *row;
  while (row = sth->fetchRow(sth)) {
    IDType id = row->getLongLongAt(row, 0);
    IDType *idP;

    if ((idP = calloc(1,sizeof(IDType))) == NULL) {
      fprintf(stderr, "Failed allocating space for a id\n");
      exit(1);
    }

    *idP = id;
    Vector_addElement(idVec, idP);
  }
  sth->finish(sth);

  Vector_setFreeFunc(idVec, free);
  return idVec;
}

/*
=head2 fetch_all_by_Transcript

  Arg[1]      : Bio::EnsEMBL::Transcript Transcript to search with
  Example     : my $ises = $isea->fetch_all_by_Transcript($transcript);
  Description : Uses the given Transcript to search for all instances of
                IntronSupportingEvidence linked to the transcript in the 
                database 
  Returntype  : ArrayRef of IntronSupportingEvidence objects
  Exceptions  : Thrown if arguments are not as stated and for DB errors 

=cut
*/
Vector *IntronSupportingEvidenceAdaptor_fetchAllByTranscript(IntronSupportingEvidenceAdaptor *isea, Transcript *transcript) {
  char qStr[1024];

  sprintf(qStr,"SELECT intron_supporting_evidence_id "
                 "FROM transcript_intron_supporting_evidence "
                "WHERE transcript_id = "IDFMTSTR, Transcript_getDbID(transcript));

  StatementHandle *sth = isea->prepare((BaseAdaptor *)isea,qStr,strlen(qStr));
  sth->execute(sth);

  Vector *idVec = Vector_new();
  ResultRow *row;
  while (row = sth->fetchRow(sth)) {
    IDType id = row->getLongLongAt(row, 0);
    IDType *idP;

    if ((idP = calloc(1,sizeof(IDType))) == NULL) {
      fprintf(stderr, "Failed allocating space for a id\n");
      exit(1);
    }

    *idP = id;
    Vector_addElement(idVec, idP);
  }
  sth->finish(sth);

  Vector *out;
  if (Vector_getNumElement(idVec) > 0) {
    out = IntronSupportingEvidenceAdaptor_fetchAllByDbIDList(isea, idVec, NULL); 
  } else {
    out = Vector_new();
  }
 
  // Free ids vector
  Vector_setFreeFunc(idVec, free);
  Vector_free(idVec);

  return out;
}

/*
=head2 fetch_flanking_exon_ids

  Arg[1]      : Bio::EnsEMBL::IntronSupportingEvidence Evidence to search with
  Arg[2]      : Bio::EnsEMBL::Transcript Transcript to search with
  Example     : my ($prev_id, $next_id) = @{$isea->fetch_flanking_exon_ids($ise, $transcript)};
  Description : Uses the given IntronSupportingEvidence and Transcript to search
                for the recorded previous and next exon database ids 
  Returntype  : ArrayRef 1 row long but with 2 columns representing previous 
                and next IDs respectivly
  Exceptions  : Thrown if arguments are not as stated and for DB errors 

=cut
*/
// Switch to passing two element array to fill in as an argument, so don't have to allocate it
IDType *IntronSupportingEvidenceAdaptor_fetchFlankingExonIds(IntronSupportingEvidenceAdaptor *isea, IntronSupportingEvidence *ise, Transcript *transcript, IDType *flanks) {
  char qStr[1024];

  sprintf(qStr,"SELECT previous_exon_id, next_exon_id "
                 "FROM transcript_intron_supporting_evidence "
                "WHERE transcript_id ="IDFMTSTR" and intron_supporting_evidence_id ="IDFMTSTR, 
          Transcript_getDbID(transcript), IntronSupportingEvidence_getDbID(ise));

  StatementHandle *sth = isea->prepare((BaseAdaptor *)isea,qStr,strlen(qStr));
  sth->execute(sth);

  if (sth->numRows(sth) == 0) {
    return NULL;
  }
  
  ResultRow *row = sth->fetchRow(sth);

  flanks[0] = row->getLongLongAt(row, 0);
  flanks[1] = row->getLongLongAt(row, 1);

  sth->finish(sth);

  return flanks;
}

Vector *IntronSupportingEvidenceAdaptor_objectsFromStatementHandle(IntronSupportingEvidenceAdaptor *isea, 
                                                                   StatementHandle *sth,
                                                                   AssemblyMapper *assMapper,
                                                                   Slice *destSlice) {
  SliceAdaptor *sa     = DBAdaptor_getSliceAdaptor(isea->dba);
  AnalysisAdaptor *aa  = DBAdaptor_getAnalysisAdaptor(isea->dba);

  Vector *features = Vector_new();
  IDHash *sliceHash = IDHash_new(IDHASH_SMALL);
  
/* Unneccesary
  my %analysis_hash;
  my %sr_name_hash;
  my %sr_cs_hash;
*/
  

  
/* Unused
  my $asm_cs;
  my $cmp_cs;
  my $asm_cs_vers;
  my $asm_cs_name;
  my $cmp_cs_vers;
  my $cmp_cs_name;
  if($mapper) {
    $asm_cs = $mapper->assembled_CoordSystem();
    $cmp_cs = $mapper->component_CoordSystem();
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
  //CoordSystem *destSliceCs;
  char *       destSliceSrName;
  IDType       destSliceSrId = 0;
  //AssemblyMapperAdaptor *asma;

  if (destSlice) {
    destSliceStart  = Slice_getStart(destSlice);
    destSliceEnd    = Slice_getEnd(destSlice);
    destSliceStrand = Slice_getStrand(destSlice);
    destSliceLength = Slice_getLength(destSlice);
    //??destSliceCs     = Slice_getCoordSystem(destSlice);
    destSliceSrName = Slice_getSeqRegionName(destSlice);
    destSliceSrId   = Slice_getSeqRegionId(destSlice);
    //??asma            = DBAdaptor_getAssemblyMapperAdaptor(ea->dba);
  }

  ResultRow *row;
  while (row = sth->fetchRow(sth)) {
    IDType id =           row->getLongLongAt(row,0);
    IDType analysisId =   row->getLongLongAt(row,1);
    IDType seqRegionId =  row->getLongLongAt(row,2);
    long seqRegionStart = row->getLongAt(row,3);
    long seqRegionEnd =   row->getLongAt(row,4);
    int seqRegionStrand = row->getIntAt(row,5);
    char *hitName =       row->getStringAt(row,6);
    double score =        row->getDoubleAt(row,7);
    char *scoreType =     row->getStringAt(row,8);
    int spliceCanonical = row->getIntAt(row,9); 

    // get the analysis object
    Analysis *analysis = AnalysisAdaptor_fetchByDbID(aa, analysisId);

/*
    // need to get the internal_seq_region, if present
    $seq_region_id = $self->get_seq_region_id_internal($seq_region_id);
    #get the slice object
    my $slice = $slice_hash{"ID:".$seq_region_id};
    if(!$slice) {
      $slice = $sa->fetch_by_seq_region_id($seq_region_id);
      $slice_hash{"ID:".$seq_region_id} = $slice;
      $sr_name_hash{$seq_region_id} = $slice->seq_region_name();
      $sr_cs_hash{$seq_region_id} = $slice->coord_system();
    }

    my $sr_name = $sr_name_hash{$seq_region_id};
    my $sr_cs   = $sr_cs_hash{$seq_region_id};
*/
    if (! IDHash_contains(sliceHash, seqRegionId)) {
      IDHash_add(sliceHash, seqRegionId, SliceAdaptor_fetchBySeqRegionId(sa, seqRegionId, POS_UNDEF, POS_UNDEF, STRAND_UNDEF));
    }
    Slice *slice = IDHash_getValue(sliceHash, seqRegionId);

    Slice *iseSlice = slice;
    
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

      
/* Was - but identical if and else so why test???
      #get a slice in the coord system we just mapped to
      if($asm_cs == $sr_cs || ($cmp_cs != $sr_cs && $asm_cs->equals($sr_cs))) {
        $slice = $slice_hash{"ID:".$seq_region_id} ||=
          $sa->fetch_by_seq_region_id($seq_region_id);
      } else {
        $slice = $slice_hash{"ID:".$seq_region_id} ||=
          $sa->fetch_by_seq_region_id($seq_region_id);
      }
*/
// Instead...
      if (! IDHash_contains(sliceHash, seqRegionId)) {
        IDHash_add(sliceHash, seqRegionId, SliceAdaptor_fetchBySeqRegionId(sa, seqRegionId, POS_UNDEF, POS_UNDEF, STRAND_UNDEF));
      }
      iseSlice = IDHash_getValue(sliceHash, seqRegionId);
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
      iseSlice = destSlice;
    }
    
    IntronSupportingEvidence *ise = IntronSupportingEvidence_new();

    IntronSupportingEvidence_setStart             (ise, seqRegionStart);
    IntronSupportingEvidence_setEnd               (ise, seqRegionEnd);
    IntronSupportingEvidence_setStrand            (ise, seqRegionStrand);
    IntronSupportingEvidence_setSlice             (ise, iseSlice);
    IntronSupportingEvidence_setAnalysis          (ise, analysis);
    IntronSupportingEvidence_setAdaptor           (ise, (BaseAdaptor *)isea);
    IntronSupportingEvidence_setDbID              (ise, id);
    IntronSupportingEvidence_setHitName           (ise, hitName);
    IntronSupportingEvidence_setScore             (ise, score);
    IntronSupportingEvidence_setScoreType         (ise, scoreType);
    IntronSupportingEvidence_setIsSpliceCanonical(ise, spliceCanonical);

    Vector_addElement(features, ise);
  }
  
  return features;
}


/*
=head2 store

  Arg[1]      : Bio::EnsEMBL::IntronSupportingEvidence Evidence to store
  Example     : $isea->store($ise);
  Description : Stores the IntronSupportingEvidence in the database. Duplicates
                are ignored.
  Returntype  : Integer The assigned database identifier
  Exceptions  : Thrown if the given object is not a IntronSupportingEvidence, 
                and for any DB exception. 

=cut
*/
// This method was a mess - tidied and rearranged
IDType IntronSupportingEvidenceAdaptor_store(IntronSupportingEvidenceAdaptor *isea, IntronSupportingEvidence *sf)  {
  if (sf == NULL) {
    fprintf(stderr,"sf is NULL in IntronSupportingEvidenceAdaptor_store\n");
    exit(1);
  }

  Class_assertType(CLASS_INTRONSUPPORTINGEVIDENCE, sf->objectType);
  
  DBAdaptor *db = isea->dba;
  
  if (! IntronSupportingEvidence_isStored(sf, db)) {
    fprintf(stderr,"ISE already stored\n");
    return IntronSupportingEvidence_getDbID(sf);
  }
  
  Analysis *analysis = IntronSupportingEvidence_getAnalysis(sf);

  if (!Analysis_isStored(analysis, db)) {
    AnalysisAdaptor_store(analysisAdaptor, analysis);
  }
  IDType analysisId = Analysis_getDbID(analysis);

  // Think the above is equivalent to this horror
  //my $analysis_id = $analysis->is_stored($db) ? $analysis->dbID() : $db->get_AnalysisAdaptor()->store($analysis);
  
/* No transfer (see GeneAdaptor for why)
  my $seq_region_id;
  ($sf, $seq_region_id) = $self->_pre_store($sf);
*/

  IDType seqRegionId = BaseFeatureAdaptor_preStore((BaseFeatureAdaptor *)isea, sf);

  char qStr[1024];
  sprintf(qStr, "INSERT IGNORE INTO intron_supporting_evidence "
                "(analysis_id, seq_region_id, seq_region_start, seq_region_end, seq_region_strand, hit_name, score, score_type, is_splice_canonical) "
                "VALUES ("IDFMTSTR","IDFMTSTR",%ld,%ld,%d,'%s',%f,'%s',%d)", 
                analysisId,
                seqRegionId,
                IntronSupportingEvidence_getSeqRegionStart(sf),
                IntronSupportingEvidence_getSeqRegionEnd(sf),
                IntronSupportingEvidence_getSeqRegionStrand(sf),
                IntronSupportingEvidence_getHitName(sf),
                IntronSupportingEvidence_getScore(sf),
                IntronSupportingEvidence_getScoreType(sf),
                IntronSupportingEvidence_getIsSpliceCanonical(sf));

  StatementHandle *sth = isea->prepare((BaseAdaptor *)isea,qStr,strlen(qStr));
  
  sth->execute(sth);
  IDType sfId = sth->getInsertId(sth);
  sth->finish(sth);

  if (!sfId) {
    sprintf(qStr,"SELECT intron_supporting_evidence_id "
                   "FROM intron_supporting_evidence "
                  "WHERE analysis_id = "IDFMTSTR
                   " AND seq_region_id = "IDFMTSTR
                   " AND seq_region_start = %ld"
                   " AND seq_region_end = %ld" 
                   " AND seq_region_strand = %d"
                   " AND hit_name = '%s'",
                analysisId,
                seqRegionId,
                IntronSupportingEvidence_getSeqRegionStart(sf),
                IntronSupportingEvidence_getSeqRegionEnd(sf),
                IntronSupportingEvidence_getSeqRegionStrand(sf),
                IntronSupportingEvidence_getHitName(sf));

    sth = isea->prepare((BaseAdaptor *)isea,qStr,strlen(qStr));
    sth->execute(sth);
    if (sth->numRows(sth) > 0) {
      ResultRow *row = sth->fetchRow(sth);
      sfId = row->getLongLongAt(row, 0);
    }
  }
  
  IntronSupportingFeature_setAdaptor(sf, isea);
  IntronSupportingFeature_setDbID(sf, sfId);

  return IntronSupportingFeature_getDbID(sf);
}

/* NIY
=head2 store_transcript_linkage

  Arg[1]      : Bio::EnsEMBL::IntronSupportingEvidence Evidence to link
  Arg[2]      : Bio::EnsEMBL::Transcript Transcript to link
  Arg[3]      : Integer an optional ID to give if the Transcript's own ID is possibly incorrect
  Example     : $isea->store_transcript_linkage($ise, $transcript);
                $isea->store_transcript_linkage($ise, $transcript, $tid);
  Description : Links a Transcript to a portion of Intron evidence
  Returntype  : None
  Exceptions  : Thrown if the given object is not a Transcript, if the 
                transcript is not stored, if the supporting evidence is not
                stored and for any DB exception. 

=cut
*/
void IntronSupportingEvidenceAdaptor_storeTranscriptLinkage(IntronSupportingEvidenceAdaptor *isea, 
                                                            IntronSupportingEvidence *sf,
                                                            Transcript *transcript,
                                                            IDType transcriptId) {
 
  if (sf == NULL || transcript == NULL) {
    fprintf(stderr,"sf or transcript is NULL in IntronSupportingEvidenceAdaptor_storeTranscriptLinkage\n");
    exit(1);
  }

  Class_assertType(CLASS_INTRONSUPPORTINGEVIDENCE, sf->objectType);
  Class_assertType(CLASS_TRANSCRIPT, transcript->objectType);
  
  if (! IntronSupportingEvidence_isStored(sf, isea->dba)) {
    fprintf(stderr,"Cannot perform the link. The IntronSupportingEvidence must be persisted first\n");
    exit(1);
  }

// Moved up so can use in sprintf
  Intron *intron = IntronSupportingEvidence_getIntron(sf, transcript);
  Exon *prevExon = Intron_getPrevExon(intron);
  Exon *nextExon = Intron_getNextExon(intron);

  if (!transcriptId) {
    transcriptId = Transcript_getDbID(transcript);
  }

  char qStr[1024];
  sprintf(qStr, "INSERT IGNORE INTO transcript_intron_supporting_evidence "
          "(transcript_id, intron_supporting_evidence_id, previous_exon_id, next_exon_id) "
          "VALUES ("IDFMTSTR","IDFMTSTR","IDFMTSTR","IDFMTSTR")", 
          transcriptId, IntronSupportingEvidence_getDbID(sf), Exon_getDbID(prevExon), Exon_getDbID(nextExon));

  StatementHandle *sth = isea->prepare((BaseAdaptor *)isea,qStr,strlen(qStr));
  
  sth->execute(sth);
  sth->finish(sth);
  
  return;
}


/* NIY
=head2 update

  Arg[1]      : Bio::EnsEMBL::IntronSupportingEvidence Evidence to update
  Example     : $isea->update($ise);
  Description : Updates all attributes of an evidence object
  Returntype  : None
  Exceptions  : Thrown if the given object is not a IntronSupportingEvidence,
                if the object is not stored and for normal DB errors

=cut
*/
/*
sub update {
  my ($self, $sf) = @_;
  assert_ref($sf, 'Bio::EnsEMBL::IntronSupportingEvidence', 'intron_supporting_evidence');
  if (! $sf->is_stored($self->db())) {
    throw "Cannot update the supporting evidence if it has not already been stored in this database";
  }
  
  my $params = [
    [$sf->analysis()->dbID(), SQL_INTEGER],
    [$sf->slice()->get_seq_region_id(), SQL_INTEGER],
    [$sf->start(), SQL_INTEGER],
    [$sf->end(), SQL_INTEGER],
    [$sf->strand(), SQL_INTEGER],
    [$sf->hit_name(), SQL_VARCHAR],
    [$sf->score(), SQL_FLOAT],
    [$sf->score_type(), SQL_VARCHAR],
    [$sf->is_splice_canonical() || 0, SQL_INTEGER],
    [$sf->dbID(), SQL_INTEGER],
  ];
  
  my $sql = <<'SQL';
UPDATE intron_supporting_evidence
SET analysis_id =?, seq_region_id =?, seq_region_start =?, 
seq_region_end =?, seq_region_strand =?, hit_name =?, score =?, score_type =?,
is_splice_canonical =?
WHERE intron_supporting_evidence_id =?
SQL
  
  $self->dbc()->sql_helper()->execute_update(-SQL => $sql, -PARAMS => $params);
  return;
}
*/


/* NIY
=head2 remove

  Arg[1]      : Bio::EnsEMBL::IntronSupportingEvidence
  Example			: $isea->remove($ise);
  Description	: Deletes the given IntronSupportingEvidence from the database. 
                This can only occur if the object has no linked transcripts
  Returntype 	: None
  Exceptions 	: Thrown if the IntronSupportingEvidence is not stored, if
                the object has linked transcripts and in the event of any
                database error

=cut

sub remove {
  my ($self, $sf) = @_;
  assert_ref($sf, 'Bio::EnsEMBL::IntronSupportingEvidence', 'intron_supporting_evidence');
  if (! $sf->is_stored($self->db())) {
    throw "Cannot delete the supporting evidence if it has not already been stored in this database";
  }
  if($sf->has_linked_transcripts()) {
    throw sprintf('Cannot delete supporting evidence %d. It still has transcripts attached', $sf->dbID());
  }
  $self->dbc()->sql_helper()->execute_update(
    -SQL => 'DELETE from intron_supporting_evidence where intron_supporting_evidence_id =?', 
    -PARAMS => [[$sf->dbID(), SQL_INTEGER]],
  );
  return;
}

=head2 remove_all_transcript_linkages

  Arg[1]      : Bio::EnsEMBL::IntronSupportingEvidence
  Example     : $isea->remove_all_transcript_linkages($ise);
  Description : Deletes the transcript links to the given IntronSupportingEvidence
  Returntype  : None
  Exceptions  : See remove_transcript_linkage

=cut

sub remove_all_transcript_linkages {
  my ($self, $sf) = @_;
  foreach my $transcript_id (@{$self->list_linked_transcript_ids($sf)}) {
    $self->_remove_transcript_linkage($sf, $transcript_id);
  }
  return;
}

=head2 remove_transcript_linkage

  Arg[1]      : Bio::EnsEMBL::IntronSupportingEvidence Evidence to unlink
  Arg[2]      : Bio::EnsEMBL::Transcript Transcript to unlink
  Example     : $isea->remove_transcript_linkages($ise, $transcript);
  Description : Deletes a transcript's link to the given IntronSupportingEvidence
  Returntype  : None
  Exceptions  : Thrown if the given object is not a Transcript, if the 
                transcript is not stored, if the supporting evidence is not
                stored and for any DB exception. 

=cut

sub remove_transcript_linkage {
  my ($self, $sf, $transcript) = @_;
  assert_ref($transcript, 'Bio::EnsEMBL::Transcript', 'transcript');
  if (! $transcript->is_stored($self->db())) {
    throw "Cannot delete the supporting evidence to transcript linkage if the transcript has not already been stored in this database";
  }
  $self->_remove_transcript_linkage($sf, $transcript->dbID());
  return;
}

sub _remove_transcript_linkage {
  my ($self, $sf, $transcript_id) = @_;
  if (! $sf->is_stored($self->db())) {
    throw "Cannot delete the supporting evidence to transcript linkage if the evidence has not already been stored in this database";
  }
  $self->dbc()->sql_helper()->execute_update(
    -SQL => 'DELETE from transcript_intron_supporting_evidence where intron_supporting_evidence_id =? and transcript_id =?', 
    -PARAMS => [[$sf->dbID(), SQL_INTEGER], [$transcript_id, SQL_INTEGER]],
  );
  return;
}
*/

