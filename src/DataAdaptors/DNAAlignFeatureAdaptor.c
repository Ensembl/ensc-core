#include "DNAAlignFeatureAdaptor.h"


char ***DNAAlignFeatureAdaptor_tableNames = {{"dna_align_feature","daf"}};

DNAAlignFeatureAdaptor *DNAAlignFeatureAdaptor_new(DBAdaptor *dba) {
  DNAAlignFeatureAdaptor *dafa;

  if ((dafa = (DNAAlignFeatureAdaptor *)calloc(1,sizeof(DNAAlignFeatureAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for DNAAlignFeatureAdaptor\n");
    return NULL;
  }
  BaseFeatureAdaptor_init((BaseFeatureAdaptor *)dafa, dba, DNAALIGNFEATURE_ADAPTOR);

  return dafa;
}

int DNAAlignFeatureAdaptor_store(BaseFeatureAdaptor *baf, Set *features) {
/*
  my ($self, @sf) = @_;

  my @tabs = $self->_tables;
  my ($tablename) = @{$tabs[0]};
  
  if( scalar(@sf) == 0 ) {
    $self->throw("Must call store with sequence features");
  }
  
  my $sth = $self->prepare("
     INSERT INTO $tablename (contig_id, contig_start, contig_end,
                             contig_strand, hit_start, hit_end,
                             hit_strand, hit_name, cigar_line,
                             analysis_id, score, evalue, perc_ident) 
     VALUES (?,?,?,?,?,?,?,?,?,?,?, ?, ?)");

  foreach my $sf ( @sf ) {
    if( !ref $sf || !$sf->isa("Bio::EnsEMBL::DnaDnaAlignFeature") ) {
      $self->throw("feature must be a Bio::EnsEMBL::DnaDnaAlignFeature," 
                    . " not a [$sf]");
    }
    
    my $contig = $sf->entire_seq();
    unless(defined $contig && $contig->isa("Bio::EnsEMBL::RawContig")) {
      $self->throw("A contig must be attached to the features to be " .
                   "stored via the attach seq method\n");
    }

    if( !defined $sf->analysis ) {
      $self->throw("Cannot store sequence features without analysis");
    }

     # will only store if object is not already stored in this database
    $self->db()->get_AnalysisAdaptor()->store( $sf->analysis() );

    $sth->execute( $contig->dbID(), $sf->start, $sf->end, $sf->strand,
                   $sf->hstart, $sf->hend, $sf->hstrand, $sf->hseqname,
                   $sf->cigar_string, $sf->analysis->dbID, $sf->score, 
                   $sf->p_value, $sf->percent_id);
    $sf->dbID($sth->{'mysql_insertid'});
  }
*/
}


char ***DNAAlignFeatureAdaptor_getTables() {
  return DNAAlignFeatureAdaptor_tableNames;
}

char *DNAAlignFeatureAdaptor_getColumns() {
  return "daf.dna_align_feature_id,"
         "daf.contig_id,"
         "daf.analysis_id,"
         "daf.contig_start,"
         "daf.contig_end,"
         "daf.contig_strand,"
         "daf.hit_start,"
         "daf.hit_end,"
         "daf.hit_name,"
         "daf.hit_strand,"
         "daf.cigar_line,"
         "daf.evalue,"
         "daf.perc_ident,"
         "daf.score";
}

Set *DNAAlignFeatureAdaptor_objectsFromStatementHandle(BaseFeatureAdaptor *baf,
                                                       StatementHandle *sth) {
 my ($self, $sth, $mapper, $slice) = @_;

  my ($dna_align_feature_id, $contig_id, $analysis_id, $contig_start, 
      $contig_end, $contig_strand, $hit_start, $hit_end, $hit_name, 
      $hit_strand, $cigar_line, $evalue, $perc_ident, $score);
  
  my $rca = $self->db()->get_RawContigAdaptor();
  my $aa = $self->db()->get_AnalysisAdaptor();
  
  my ($analysis, $contig);
  my @features;

  my %a_hash;

  my ($row, $row_cache);

  $row_cache = $sth->fetchall_arrayref();

  if($slice) {
    my ($chr, $start, $end, $strand);
    my $slice_start  = $slice->chr_start();
    my $slice_end    = $slice->chr_end();
    my $slice_strand = $slice->strand();
    my $slice_name   = $slice->name();

    my ($feat_start, $feat_end, $feat_strand);

    while($row = shift @$row_cache) {
      ($dna_align_feature_id, $contig_id, $analysis_id, $contig_start, 
       $contig_end, $contig_strand, $hit_start, $hit_end, $hit_name, 
       $hit_strand, $cigar_line, $evalue, $perc_ident, $score) = @$row;

      #convert contig coordinates to assembly coordinates
      ($chr, $start, $end, $strand) = 
        $mapper->fast_to_assembly($contig_id, $contig_start, 
                                  $contig_end, $contig_strand);
      
      #if mapped to gap, skip
      next unless(defined $start);

      #if mapped outside slice region, skip
      next if ($start > $slice_end) || ($end < $slice_start); 

      #convert assembly coordinates to slice coordinates
      if($slice_strand == -1) {
        $feat_start  = $slice_end - $end + 1;
        $feat_end    = $slice_end - $start + 1;
        $feat_strand = $strand * -1;
      } else {
        $feat_start  = $start - $slice_start + 1;
        $feat_end    = $end   - $slice_start + 1;
        $feat_strand = $strand;
      }

      $analysis = $a_hash{$analysis_id} ||= $aa->fetch_by_dbID($analysis_id);

      push @features, Bio::EnsEMBL::DnaDnaAlignFeature->new_fast(
                    {'_gsf_tag_hash'  =>  {},
                     '_gsf_sub_array' =>  [],
                     '_parse_h'       =>  {},
                     '_analysis'      =>  $analysis,
                     '_gsf_start'     =>  $feat_start,
                     '_gsf_end'       =>  $feat_end,
                     '_gsf_strand'    =>  $feat_strand,
                     '_gsf_score'     =>  $score,
                     '_seqname'       =>  $slice_name,
                     '_percent_id'    =>  $perc_ident,
                     '_p_value'       =>  $evalue,
                     '_hstart'        =>  $hit_start,
                     '_hend'          =>  $hit_end,
                     '_hstrand'       =>  $hit_strand,
                     '_hseqname'      =>  $hit_name,
                     '_gsf_seq'       =>  $slice,
                     '_cigar_string'  =>  $cigar_line,
                     '_id'            =>  $hit_name,
                     '_database_id'   =>  $dna_align_feature_id});
    }
  } else {
    my %c_hash;
    while($row = shift @$row_cache) {
      ($dna_align_feature_id, $contig_id, $analysis_id, $contig_start, 
       $contig_end, $contig_strand, $hit_start, $hit_end, $hit_name, 
       $hit_strand, $cigar_line, $evalue, $perc_ident, $score) = @$row;
      
      $analysis = $a_hash{$analysis_id} ||= $aa->fetch_by_dbID($analysis_id);
      $contig   = $c_hash{$contig_id}   ||= $rca->fetch_by_dbID($contig_id);

      #use a very fast (hack) constructor since we may be creating over 10000
      #features at a time and normal object construction is too slow.
      push @features, Bio::EnsEMBL::DnaDnaAlignFeature->new_fast(
                    {'_gsf_tag_hash'  =>  {},
                     '_gsf_sub_array' =>  [],
                     '_parse_h'       =>  {},
                     '_analysis'      =>  $analysis,
                     '_gsf_start'     =>  $contig_start,
                     '_gsf_end'       =>  $contig_end,
                     '_gsf_strand'    =>  $contig_strand,
                     '_gsf_score'     =>  $score,
                     '_seqname'       =>  $contig->name,
                     '_percent_id'    =>  $perc_ident,
                     '_p_value'       =>  $evalue,
                     '_hstart'        =>  $hit_start,
                     '_hend'          =>  $hit_end,
                     '_hstrand'       =>  $hit_strand,
                     '_hseqname'      =>  $hit_name,
                     '_gsf_seq'       =>  $contig,
                     '_cigar_string'  =>  $cigar_line,
                     '_id'            =>  $hit_name,
                     '_database_id'   =>  $dna_align_feature_id}); 

    }
    
  }
  
  return \@features;
}
