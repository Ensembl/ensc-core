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


=head2 store

  Arg [1]    : list of Bio::EnsEMBL::SimpleFeatures @sf
               the simple features to store in the database
  Example    : $simple_feature_adaptor->store(1234, @simple_feats);
  Description: Stores a list of simple feature objects in the database
  Returntype : none
  Exceptions : thrown if @sf is not defined, if any of the features do not
               have an attached contig object, 
               or if any elements of @sf are not Bio::EnsEMBL::SeqFeatures 
  Caller     : general

=cut

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


=head2 _tablename

  Arg [1]    : none
  Example    : none
  Description: PROTECTED implementation of superclass abstract method
               returns the names, aliases of the tables to use for queries
  Returntype : list of listrefs of strings
  Exceptions : none
  Caller     : internal

=cut

char ***SimpleFeatureAdaptor_getTables() {
  return SimpleFeatureAdaptor_tableNames;
}


=head2 _columns

  Arg [1]    : none
  Example    : none
  Description: PROTECTED implementation of superclass abstract method
               returns a list of columns to use for queries
  Returntype : list of strings
  Exceptions : none
  Caller     : internal

=cut

char *SimpleFeatureAdaptor_getColumns() {
  return "sf.simple_feature_id " 
	 "sf.contig_id "
         "sf.contig_start "
         "sf.contig_end "
         "sf.contig_strand "
	 "sf.display_label "
         "sf.analysis_id "
         "sf.score";
}


=head2 _objs_from_sth

  Arg [1]    : hash reference $hashref
  Example    : none
  Description: PROTECTED implementation of superclass abstract method.
               creates SimpleFeatures from an executed DBI statement handle.
  Returntype : list reference to Bio::EnsEMBL::SimpleFeature objects
  Exceptions : none
  Caller     : internal

=cut

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
