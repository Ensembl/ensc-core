#include "RepeatFeatureAdaptor.h"

#include "IDHash.h"
#include "RepeatConsensusAdaptor.h"
#include "AnalysisAdaptor.h"
#include "RawContigAdaptor.h"

#include "Repeat.h"
#include "RepeatConsensus.h"

NameTableType RepeatFeatureAdaptor_tableNames = {{"repeat_feature","r"},
                                                 {"repeat_consensus","rc"},
                                                 {NULL,NULL}};

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

NameTableType *RepeatFeatureAdaptor_getTables() {
  return &RepeatFeatureAdaptor_tableNames;
}

char *RepeatFeatureAdaptor_getColumns() {
  return "r.repeat_feature_id,"
         "r.contig_id,"
         "r.analysis_id,"
         "r.contig_start,"
         "r.contig_end,"
         "r.contig_strand,"
         "r.repeat_consensus_id,"
         "r.repeat_start,"
         "r.repeat_end,"
         "r.score,"
         "rc.repeat_name,"
         "rc.repeat_class,"
         "rc.repeat_consensus";
}

Set *RepeatFeatureAdaptor_objectsFromStatementHandle(BaseFeatureAdaptor *baf,
                                                     StatementHandle *sth,
                                                     AssemblyMapper *mapper,
                                                     Slice *slice) {
  AnalysisAdaptor *aa;
  RawContigAdaptor *rca;
  RepeatConsensusAdaptor *rpca;
  Set *features;
  ResultRow *row;
  int i;
  IDHash *rcHash;

  aa = DBAdaptor_getAnalysisAdaptor(bfa->dba);
  rca = DBAdaptor_getRawContigAdaptor(bfa->dba);
  rpca = DBAdaptor_getRepeatConsensusAdaptor(bfa->dba);


  features = Set_new();
  rcHash = IDHash_new(IDHASH_SMALL);


  my $rc;
  my $contig;
  my $analysis;

  while (row = sth->fetchRow(sth)) {
    RepeatFeature *rf;
    Analysis  *analysis = AnalysisAdaptor_fetchByDbID(aa, row->getLongLongAt(row,2));
    RawContig *contig = RawContigAdaptor_fetchByDbID(rca, row->getLongLongAt(row,1));
    int64 repeatConsensusId = row->getLongLongAt(row,6);
    RepeatConsensus *rc;

  $sth->bind_columns( \$repeat_feature_id, \$contig_id, \$analysis_id, \$contig_start, 
                      \$contig_end, \$contig_strand, \$repeat_consensus_id, 
                      \$repeat_start,\$repeat_end, \$score, 
                      \$repeat_name, \$repeat_class,
                      \$repeat_consensus );

    //create a repeat consensus object
    if (!IDHash_constains(rcHash, repeatConsensusId)) {
      rc = RepeatConsensus_new();
      RepeatConsensus_setDbID(rc, repeatConsensusId);
      RepeatConsensus_setName(rc, row->getStringAt(10));
      RepeatConsensus_setRepeatClass(rc, row->getStringAt(11));
      RepeatConsensus_setConsensus(rc, row->getStringAt(12));
      RepeatConsensus_setAdaptor(rpca);

      IDHash_add(rcHash,repeatConsensusId);
    } else {
      rc = IDHash_get(rcHash,repeatConsensusId);
    }
    
    //create the new repeat feature
    push @features, Bio::EnsEMBL::RepeatFeature->new_fast(
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
    Set_addElement(features,rf);
  }

  IDHash_free(rcHash);

  return features;
}

char *RepeatFeatureAdaptor_defaultWhereClause() {
  return "r.repeat_consensus_id = rc.repeat_consensus_id";
}

