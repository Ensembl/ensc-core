#include "PredictionTranscriptAdaptor.h"


NameTableType PredictionTranscriptAdaptor_tableNames = {{"prediction_transcript","p"},
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
  pta->finalClause = PredictionTranscriptAdaptor_finalClause;

  return pta;
}

int PredictionTranscriptAdaptor_store(BaseFeatureAdaptor *bfa, Set *features) {
/*
  my ( $self, @pre_transcripts ) = @_;

  my $exon_sql = q{
      INSERT INTO prediction_transcript ( prediction_transcript_id, exon_rank,
                                          contig_id, contig_start, contig_end,
                                          contig_strand, start_phase, score,
                                          p_value, analysis_id, exon_count )
        VALUES ( ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ? )
      };

  my $exonst = $self->prepare($exon_sql);

  foreach my $pre_trans (@pre_transcripts) {
    if( ! $pre_trans->isa('Bio::EnsEMBL::PredictionTranscript') ) {
      $self->throw("$pre_trans is not a EnsEMBL PredictionTranscript "
                   . "- not dumping!");
    }

    if( $pre_trans->dbID && $pre_trans->adaptor == $self ) {
      $self->warn("Already stored");
    }

    my $exonId = undef;
    my $exons = $pre_trans->get_all_Exons();
    my $dbID = undef;
    my $rank = 1;

    for my $exon ( @$exons ) {
      if( ! defined $exon ) { $rank++; next; }

      my $contig_id = $exon->contig->dbID();
      my $contig_start = $exon->start();
      my $contig_end = $exon->end();
      my $contig_strand = $exon->strand();

      my $start_phase = $exon->phase();
      my $end_phase = $exon->end_phase();

      # this is only in PredictionExon
      my $score = $exon->score();
      my $p_value = $exon->p_value();

      my $analysis = $pre_trans->analysis->dbID;

      if( $rank == 1 ) {
        $exonst->execute( undef, 1, $contig_id, $contig_start,
                          $contig_end, $contig_strand,
                          $start_phase, $score, $p_value, $analysis,
                          scalar( @{$exons} ));
        $dbID = $exonst->{'mysql_insertid'};
      } else {
        $exonst->execute( $dbID, $rank, $contig_id, $contig_start,
                          $contig_end, $contig_strand,
                          $start_phase, $score, $p_value, $analysis,
                         scalar( @{$exons} ) );
      }
      $rank++;
    }

    $pre_trans->dbID( $dbID );
    $pre_trans->adaptor( $self );
  }

  $exonst->finish;
*/
}


NameTableType *PredictionTranscriptAdaptor_getTables() {
  return &PredictionTranscriptAdaptor_tableNames;
}

char *PredictionTranscriptAdaptor_getColumns() {
  return "p.prediction_transcript_id,"
         "p.contig_id,"
         "p.contig_start,"
         "p.contig_end,"
         "p.contig_strand,"
         "p.start_phase,"
         "p.exon_rank,"
         "p.score,"
         "p.p_value,"
         "p.analysis_id,"
         "p.exon_count";
}

Set *PredictionTranscriptAdaptor_objectsFromStatementHandle(BaseFeatureAdaptor *bfa,
                                                       StatementHandle *sth,
                                                       AssemblyMapper *mapper,
                                                       Slice *slice) {
  my ($self, $sth, $mapper, $slice) = @_;

  my @out = ();

  my ($prediction_transcript_id,
      $contig_id, $contig_start, $contig_end, $contig_strand,
      $start_phase, $exon_rank, $score, $p_value, $analysis_id,
      $exon_count );

  $sth->bind_columns(\$prediction_transcript_id,
                    \$contig_id, \$contig_start, \$contig_end, \$contig_strand,
                    \$start_phase, \$exon_rank, \$score, \$p_value,
                    \$analysis_id,\$exon_count);

  my $rca = $self->db->get_RawContigAdaptor;
  my $aa  = $self->db->get_AnalysisAdaptor;

  my ($analysis, $contig, $pre_trans, $ptid, $on_slice_flag, $last_end,
      $chr, $start, $end, $strand,
      $slice_start, $slice_end, $slice_strand,
      $exon, $exon_start, $exon_end, $exon_strand,
      $stable_start, $stable_end, $stable_ctg,
      $transcript_slice_start, $transcript_slice_end );
  my (%analysis_hash, %contig_hash);

  if($slice) {
    $slice_start  = $slice->chr_start;
    $slice_end    = $slice->chr_end;
    $slice_strand = $slice->strand;
  }

  $on_slice_flag = 0;


  while($sth->fetch) {
    #create a new transcript for each new prediction transcript id
    unless(defined $pre_trans && $ptid == $prediction_transcript_id) {
      $pre_trans = Bio::EnsEMBL::PredictionTranscript->new;

      $ptid = $prediction_transcript_id;
      $pre_trans->dbID($ptid);

      unless($analysis = $analysis_hash{$analysis_id}) {
        $analysis = $aa->fetch_by_dbID($analysis_id);
        $analysis_hash{$analysis_id} = $analysis;
      }

      $pre_trans->analysis($analysis);
      $pre_trans->set_exon_count($exon_count);

      if(@out) {
        #throw away last pt if no exons or introns were on the slice
        if($slice && ( $transcript_slice_end < 1 ||
                       $transcript_slice_start > $slice->length() )) {
          pop @out;
        } else {
          #set the stable_id of the previous prediction
          $out[$#out]->stable_id("$stable_ctg.$stable_start.$stable_end");
        }
      }

      push( @out, $pre_trans );

      #reset values used for last predtrans
      $stable_start = -1;
      $stable_end   = -1;
      $stable_ctg = '';

      $transcript_slice_end = undef;
      $transcript_slice_start = undef;
    }

    #recalculate stable id values
    if($stable_start == -1 || $contig_start < $stable_start) {
      $stable_start = $contig_start;
    }
    if($contig_end > $stable_end) {
      $stable_end = $contig_end;
    }
    unless($contig = $contig_hash{$contig_id}) {
      $contig = $rca->fetch_by_dbID($contig_id);
      $contig_hash{$contig_id} = $contig;
    }
    $stable_ctg = $contig->name;

    if($slice) {
      #a slice was passed in so we want slice coords

      #convert contig coords to assembly coords
      ($chr, $start, $end, $strand) =
        $mapper->fast_to_assembly($contig_id, $contig_start,
                                  $contig_end, $contig_strand);

      #if mapped to gap skip
      next unless(defined $start);


      #convert to slice coordinates
      if($slice_strand == -1) {
        $exon_start  = $slice_end - $end   + 1;
        $exon_end    = $slice_end - $start + 1;
        $exon_strand = $strand * -1;
      } else {
        $exon_start  = $start - $slice_start + 1;
        $exon_end    = $end   - $slice_start   + 1;
        $exon_strand = $strand;
      }  

      if( !defined $transcript_slice_start ||
          $transcript_slice_start > $exon_start ) {
        $transcript_slice_start = $exon_start;
      }

      if( ! defined $transcript_slice_end ||
          $transcript_slice_end < $exon_end ) {
        $transcript_slice_end = $exon_end;
      }
      #use slice as the contig instead of the raw contig
      $contig = $slice;
    } else {
      #we just want plain old contig coords
      $exon_start =  $contig_start;
      $exon_end   =  $contig_end;
      $exon_strand = $contig_strand;
    }

    #create an exon and add it to the prediction transcript
    $exon = Bio::EnsEMBL::Exon->new_fast($contig,
                                         $exon_start,
                                         $exon_end,
                                         $exon_strand);
    $exon->phase( $start_phase );
    $exon->end_phase( ($exon_end - $exon_start + 1 + $start_phase) % 3 );
    $exon->score( $score );
    $exon->p_value( $p_value );

    $pre_trans->add_Exon($exon, $exon_rank);
  }

  #throw away last  pred_transcript if it had no exons overlapping the slice
  if(@out) {
    if($slice && ( $transcript_slice_end < 1 ||
                   $transcript_slice_start > $slice->length() )) {
      pop @out;
    } else {
      #set the stable id of the last prediction transcript
      $out[$#out]->stable_id("$stable_ctg.$stable_start.$stable_end");
    }
  }

  return \@out;
}

char *PredictionTranscriptAdaptor_finalClause(void) {
  return  "order by p.prediction_transcript_id, p.exon_rank";
}

