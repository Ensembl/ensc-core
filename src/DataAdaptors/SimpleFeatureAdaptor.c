#include "SimpleFeatureAdaptor.h"


char ***SimpleFeatureAdaptor_tableNames = {{"simple_feature","sf"}};

SimpleFeatureAdaptor *SimpleFeatureAdaptor_new(DBAdaptor *dba) {
  SimpleFeatureAdaptor *sfa;

  if ((sfa = (SimpleFeatureAdaptor *)calloc(1,sizeof(SimpleFeatureAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for SimpleFeatureAdaptor\n");
    return NULL;
  }
  BaseFeatureAdaptor_init((BaseFeatureAdaptor *)sfa, dba, SIMPLEFEATURE_ADAPTOR);

  return sfa;
}

int SimpleFeatureAdaptor_store(BaseFeatureAdaptor *baf, Set *features) {
/*
  my ($self,@sf) = @_;
  
  if( scalar(@sf) == 0 ) {
    $self->throw("Must call store with list of sequence features");
  }
  
  my $sth = 
    $self->prepare("INSERT INTO simple_feature (contig_id, contig_start,
                                                contig_end, contig_strand,
                                                display_label, analysis_id,
                                                score) 
                    VALUES (?,?,?,?,?,?,?)");

  foreach my $sf ( @sf ) {
    if( !ref $sf || !$sf->isa("Bio::EnsEMBL::SimpleFeature") ) {
      $self->throw("Simple feature must be an Ensembl SimpleFeature, " .
		   "not a [$sf]");
    }
    
    if( !defined $sf->analysis ) {
      $self->throw("Cannot store sequence features without analysis");
    }
    if( !defined $sf->analysis->dbID ) {
      $self->throw("I think we should always have an analysis object " .
		   "which has originated from the database. No dbID, " .
		   "not putting in!");
    }
    
    my $contig = $sf->entire_seq();
    unless(defined $contig && $contig->isa("Bio::EnsEMBL::RawContig")) {
      $self->throw("Cannot store feature without a Contig object attached via "
		   . "attach_seq\n");
    }

    $sth->execute($contig->dbID(), $sf->start, $sf->end, $sf->strand,
		  $sf->display_label, $sf->analysis->dbID, $sf->score);
  } 
*/
}

char ***SimpleFeatureAdaptor_getTables() {
  return SimpleFeatureAdaptor_tableNames;
}

char *SimpleFeatureAdaptor_getColumns() {
  return "sf.simple_feature_id," 
	 "sf.contig_id,"
         "sf.contig_start,"
         "sf.contig_end,"
         "sf.contig_strand,"
	 "sf.display_label,"
         "sf.analysis_id,"
         "sf.score";
}

Set *SimpleFeatureAdaptor_objectsFromStatementHandle(BaseFeatureAdaptor *baf,
                                                     StatementHandle *sth) {
  my ($self, $sth) = @_;

  my $aa = $self->db()->get_AnalysisAdaptor();  
  my $rca = $self->db()->get_RawContigAdaptor();

  my @features = ();
  
  my $hashref;
  while($hashref = $sth->fetchrow_hashref()) {
    my $contig = $rca->fetch_by_dbID($hashref->{'contig_id'});
    my $analysis = $aa->fetch_by_dbID($hashref->{'analysis_id'});

    my $out = Bio::EnsEMBL::SimpleFeature->new();
    $out->start($hashref->{'contig_start'});
    $out->end($hashref->{'contig_end'});
    $out->strand($hashref->{'contig_strand'});
    $out->analysis($analysis);
    $out->display_label($hashref->{'display_label'});
    $out->attach_seq($contig); 

    if($hashref->{'score'}) {
      $out->score($hashref->{'score'});
    }
    
    $out->dbID($hashref->{'simple_feature_id'});

    push @features, $out;
  }

  return \@features;
}
