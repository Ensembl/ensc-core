#include "RepeatFeatureAdaptor.h"


char ***RepeatFeatureAdaptor_tableNames = {{"repeat_feature","r"},
                                           {"repeat_consensus","rc"}};

RepeatFeatureAdaptor *RepeatFeatureAdaptor_new(DBAdaptor *dba) {
  RepeatFeatureAdaptor *rfa;

  if ((rfa = (RepeatFeatureAdaptor *)calloc(1,sizeof(RepeatFeatureAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for RepeatFeatureAdaptor\n");
    return NULL;
  }
  BaseFeatureAdaptor_init((BaseFeatureAdaptor *)rfa, dba, REPEATFEATURE_ADAPTOR);

  rfa->getTables = RepeatFeatureAdaptor_getTables;
  rfa->getColumns = RepeatFeatureAdaptor_getColumns;
  rfa->store = RepeatFeatureAdaptor_store;
  rfa->objectsFromStatementHandle = RepeatFeatureAdaptor_objectsFromStatementHandle;
  rfa->defaultWhereClause = RepeatFeatureAdaptor_defaultWhereClause;

  return rfa;
}

int RepeatFeatureAdaptor_store(BaseFeatureAdaptor *baf, Set *features) {
/*
  my( $self, @repeats ) = @_;
  
  my $rca = $self->db->get_RepeatConsensusAdaptor;
  my ($cons, $db_id);

  my $sth = $self->prepare(qq{
    INSERT into repeat_feature( repeat_feature_id
                        , contig_id
                        , contig_start
                        , contig_end
                        , contig_strand
                        , repeat_consensus_id
                        , repeat_start
                        , repeat_end
                        , score
                        , analysis_id )
      VALUES(NULL, ?,?,?,?,?,?,?,?,?)
    });

  foreach my $rf (@repeats) {
    $self->throw("Must have a RepeatConsensus attached")
      unless defined ($cons = $rf->repeat_consensus);
      
    # for tandem repeats - simply store consensus and repeat
    # one pair per hit. don't need to check consensi stored
    # already. consensus has name and class set to 'trf'

    if ($cons->repeat_class eq 'trf') {

      # Look for matches already stored
      my @match = @{$rca->fetch_by_class_seq('trf', $cons->repeat_consensus)}; 
      if (@match) {
      $cons->dbID($match[0]->dbID());
      }
      else {
      $rca->store($cons);
      }

    } elsif ($cons->repeat_class eq 'Simple_repeat') {

      my $rcon = $cons->name;
      $rcon =~ s/\((\S+)\)n/$1/;   # get repeat element
      $cons->repeat_consensus($rcon);

      # Look for matches already stored
      my $match = $rca->fetch_by_name_class($cons->name, 'Simple_repeat'); 
      if ($match) {
      $cons->dbID($match->dbID());
      }
      else {
      $rca->store($cons);
      }
    } else {

      # for other repeats - need to see if a consensus is stored already
      unless ($cons->dbID) {
    my $match = ($rca->fetch_by_name($cons->name));

    if($match) {
      #set the consensus dbID to be the same as the database one
      $cons->dbID($match->dbID());
    } else {
      # if we don't match a consensus already stored create a fake one 
      # and set consensus to 'N' as null seq not allowed
      # FIXME: not happy with this, but ho hum ...
      $self->warn("Can't find " . $cons->name . "\n");
      $cons->repeat_consensus("N");
      $rca->store($cons);
    }

    #if (@match > 1) {
      #multiple consensi were matched
    #  $self->warn(@match . " consensi for " . $cons->name . "\n");
    #}
      }
    }
    
    my $contig = $rf->entire_seq();

    unless(defined $contig && $contig->isa("Bio::EnsEMBL::RawContig")) {
      $self->throw("RepeatFeature cannot be stored without a contig " .
           "attached via the attach_seq method");
    } unless($contig->dbID()) {
      $self->throw("RepeatFeature cannot be stored because attached contig " .
           "does not have a dbID");
    }
    
    $sth->execute(
          $contig->dbID(),
          $rf->start,
          $rf->end,
          $rf->strand,
          $rf->repeat_consensus->dbID(),
          $rf->hstart,
          $rf->hend,
          $rf->score,
          $rf->analysis->dbID,
         );

    my $db_id = $sth->{'mysql_insertid'}
    or $self->throw("Didn't get an insertid from the INSERT statement");
    $rf->dbID($db_id);
  }
*/
}

char ***RepeatFeatureAdaptor_getTables() {
  return RepeatFeatureAdaptor_tableNames;
}

char *RepeatFeatureAdaptor_getColumns() {
  return "r.repeat_feature_id,"
         "r.contig_id,"
         "r.contig_start,"
         "r.contig_end,"
         "r.contig_strand,"
         "r.repeat_consensus_id,"
         "r.repeat_start,"
         "r.repeat_end,"
         "r.analysis_id,"
         "r.score,"
         "rc.repeat_name,"
         "rc.repeat_class,"
         "rc.repeat_consensus";
}

Set *RepeatFeatureAdaptor_objectsFromStatementHandle(BaseFeatureAdaptor *baf,
                                                     StatementHandle *sth) {
sub _objs_from_sth {
  my ($self, $sth) = @_;

  my $rca = $self->db()->get_RepeatConsensusAdaptor();
  my $ca = $self->db()->get_RawContigAdaptor();
  my $aa = $self->db->get_AnalysisAdaptor();

  my @features;
  my %rc_hash;
  my %analysis_hash;
  my %contig_hash;

  my($repeat_feature_id, $contig_id, $contig_start, $contig_end, 
     $contig_strand, $repeat_consensus_id, $repeat_start, $repeat_end,
     $analysis_id, $score, $repeat_name, $repeat_class,
     $repeat_consensus);
  
  $sth->bind_columns( \$repeat_feature_id, \$contig_id, \$contig_start, 
                      \$contig_end, \$contig_strand, \$repeat_consensus_id, 
                      \$repeat_start,\$repeat_end, \$analysis_id, \$score, 
                      \$repeat_name, \$repeat_class,
                      \$repeat_consensus );

  my $rc;
  my $contig;
  my $analysis;

  while($sth->fetch()) {
    #create a repeat consensus object
    unless($rc = $rc_hash{$repeat_consensus_id}) {
      $rc = new Bio::EnsEMBL::RepeatConsensus;
      $rc->dbID($repeat_consensus_id);
      $rc->repeat_class($repeat_class);
      $rc->name($repeat_name);
      $rc->repeat_consensus($repeat_consensus);
      $rc->adaptor($rca);

      $rc_hash{$repeat_consensus_id} = $rc;
    }
    
    unless($analysis = $analysis_hash{$analysis_id}) {
      $analysis = $aa->fetch_by_dbID($analysis_id);
      $analysis_hash{$analysis_id} = $analysis;
    }

    unless($contig = $contig_hash{$contig_id}) {
      $contig = $ca->fetch_by_dbID($contig_id);
      $contig_hash{$contig_id} = $contig;
    }

    #create the new repeat feature
    push @features, Bio::EnsEMBL::RepeatFeature->new_fast(
                    { '_gsf_tag_hash'  =>  {},
                      '_gsf_sub_array' =>  [],
                              '_parse_h'       =>  {},
                              '_analysis'      =>  $analysis,
                              '_gsf_start'         =>  $contig_start,
                              '_gsf_end'           =>  $contig_end,
                              '_gsf_strand'        =>  $contig_strand,
                              '_gsf_score'         =>  $score,
                              '_hstart'        =>  $repeat_start,
                              '_hend'          =>  $repeat_end,
                              '_repeat_consensus' => $rc,
                      '_adaptor'       =>  $self,
                      '_gsf_seq'       =>  $contig,
                              '_db_id'         =>  $repeat_feature_id } );
  }

  return \@features;
}

char *RepeatFeatureAdaptor_defaultWhereClause() {
  return "r.repeat_consensus_id = rc.repeat_consensus_id";
}

