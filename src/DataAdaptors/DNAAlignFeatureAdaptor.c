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
}


char ***DNAAlignFeatureAdaptor_getTables() {
  return DNAAlignFeatureAdaptor_tableNames;
}

char *DNAAlignFeatureAdaptor_getColumns() {
}

Set *DNAAlignFeatureAdaptor_objectsFromStatementHandle(BaseFeatureAdaptor *baf,
                                                       StatementHandle *sth) {
}
