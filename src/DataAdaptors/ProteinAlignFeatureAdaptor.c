#include "ProteinAlignFeatureAdaptor.h"

#include "RawContigAdaptor.h"
#include "AnalysisAdaptor.h"

NameTableType ProteinAlignFeatureAdaptor_tableNames = {{"protein_align_feature","paf"},{NULL,NULL}};

ProteinAlignFeatureAdaptor *ProteinAlignFeatureAdaptor_new(DBAdaptor *dba) {
  ProteinAlignFeatureAdaptor *pafa;

  if ((pafa = (ProteinAlignFeatureAdaptor *)calloc(1,sizeof(ProteinAlignFeatureAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for ProteinAlignFeatureAdaptor\n");
    return NULL;
  }
  BaseFeatureAdaptor_init((BaseFeatureAdaptor *)pafa, dba, PROTEINALIGNFEATURE_ADAPTOR);

  pafa->getTables = ProteinAlignFeatureAdaptor_getTables;
  pafa->getColumns = ProteinAlignFeatureAdaptor_getColumns;
  pafa->store = ProteinAlignFeatureAdaptor_store;
  pafa->objectsFromStatementHandle = ProteinAlignFeatureAdaptor_objectsFromStatementHandle;

  return pafa;
}

int ProteinAlignFeatureAdaptor_store(BaseFeatureAdaptor *baf, Set *features) {
/*
   my ($self ,@sf) = @_;

   if( scalar(@sf) == 0 ) {
       $self->throw("Must call store with list of sequence features");
   }

   my $sth = $self->prepare(
        "INSERT INTO protein_align_feature(contig_id, contig_start, contig_end,
                                           contig_strand, hit_start, hit_end,
                                           hit_name, cigar_line, analysis_id,
                                           score, evalue, perc_ident) 
         VALUES (?,?,?,?,?,?,?,?,?,?, ?, ?)");

   foreach my $sf ( @sf ) {
     if( !ref $sf || !$sf->isa("Bio::EnsEMBL::DnaPepAlignFeature") ) {
       $self->throw("Feature must be a Bio::EnsEMBL::DnaPepAlignFeature, " .
                    "not a [$sf]");
     }
     
     if( !defined $sf->analysis ) {
       $self->throw("Cannot store sequence features without analysis");
     }

     # will only store if object is not already stored in this database
     $self->db()->get_AnalysisAdaptor()->store( $sf->analysis() );

     my $contig = $sf->entire_seq();
     #print STDERR $contig."\n";
     unless(defined $contig && $contig->isa("Bio::EnsEMBL::RawContig")) { 
       $self->throw("Cannot store feature without Contig attached via " .
                    "attach_seq\n");
     }   
     
     $sth->execute($contig->dbID(), $sf->start, $sf->end,
                   $sf->strand, $sf->hstart, $sf->hend,
                   $sf->hseqname, $sf->cigar_string, $sf->analysis->dbID,
                   $sf->score, $sf->p_value, $sf->percent_id);
     
     $sf->dbID($sth->{'mysql_insertid'});
   }
*/

}


NameTableType *ProteinAlignFeatureAdaptor_getTables() {
  return &ProteinAlignFeatureAdaptor_tableNames;
}

char *ProteinAlignFeatureAdaptor_getColumns() {

  return "paf.protein_align_feature_id,"
         "paf.contig_id,"
         "paf.contig_start,"
         "paf.contig_end,"
         "paf.analysis_id,"
         "paf.contig_strand,"
         "paf.hit_start,"
         "paf.hit_end,"
         "paf.hit_name,"
         "paf.cigar_line,"
         "paf.evalue,"
         "paf.perc_ident,"
         "paf.score";

}

Set *ProteinAlignFeatureAdaptor_objectsFromStatementHandle(BaseFeatureAdaptor *baf,
                                                           StatementHandle *sth) {
  my ($self, $sth, $mapper, $slice) = @_;

  my ($protein_align_feature_id, $contig_id, $contig_start, $contig_end,
      $analysis_id, $contig_strand, $hit_start, $hit_end, $hit_name, 
      $cigar_line, $evalue, $perc_ident, $score);
  
  my $rca = $self->db()->get_RawContigAdaptor();
  my $aa = $self->db()->get_AnalysisAdaptor();
  my @features;

  my($analysis, $contig);

  my %a_hash;
  my %c_hash;

  my($row_cache, $row);

  $row_cache = $sth->fetchall_arrayref();

  if($slice) {
    my ($chr, $start, $end, $strand);
    my $slice_start   = $slice->chr_start();
    my $slice_end     = $slice->chr_end();
    my $slice_name    = $slice->name();
    my $slice_strand = $slice->strand();

    my($feat_start, $feat_end, $feat_strand);

    while($row = shift @$row_cache) {
      ($protein_align_feature_id, $contig_id, $contig_start, $contig_end,
      $analysis_id, $contig_strand, $hit_start, $hit_end, $hit_name, 
      $cigar_line, $evalue, $perc_ident, $score) = @$row;

      $analysis = $a_hash{$analysis_id} ||= $aa->fetch_by_dbID($analysis_id);
      ($chr, $start, $end, $strand) = 
        $mapper->fast_to_assembly($contig_id, $contig_start, 
                                  $contig_end, $contig_strand);
      
      #skip if feature maps to gap
      next unless(defined $start);

      #skip if feature outside of slice area
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

      #use a very fast (hack) constructor - normal object construction is too
      #slow for the number of features we are potentially dealing with
      push @features, Bio::EnsEMBL::DnaPepAlignFeature->new_fast(
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
                 '_hstrand'       =>  1, #strand is always one for pep hits
                 '_hseqname'      =>  $hit_name,
                 '_gsf_seq'       =>  $slice,
                 '_cigar_string'  =>  $cigar_line,
                 '_id'            =>  $hit_name,
                 '_database_id'   =>  $protein_align_feature_id});
    }    

  } else {
    while($row = shift @$row_cache) {
      ($protein_align_feature_id, $contig_id, $contig_start, $contig_end,
      $analysis_id, $contig_strand, $hit_start, $hit_end, $hit_name, 
      $cigar_line, $evalue, $perc_ident, $score) = @$row;

      $analysis = $a_hash{$analysis_id} ||= $aa->fetch_by_dbID($analysis_id);
      $contig   = $c_hash{$contig_id}   ||= $rca->fetch_by_dbID($contig_id);
      
      #use a very fast (hack) constructor - normal object construction is too
      #slow for the number of features we are potentially dealing with
      push @features, Bio::EnsEMBL::DnaPepAlignFeature->new_fast(
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
                 '_hstrand'       =>  1,  #strand is always 1 for pep hits
                 '_hseqname'      =>  $hit_name,
                 '_gsf_seq'       =>  $contig,
                 '_cigar_string'  =>  $cigar_line,
                 '_id'            =>  $hit_name,
                 '_database_id'   =>  $protein_align_feature_id});
    }
  }


  return \@features;
}
