#include "RepeatConsensusAdaptor.h"
#include "BaseAdaptor.h"

#include "StatementHandle.h"
#include "ResultRow.h"


RepeatConsensusAdaptor *RepeatConsensusAdaptor_new(DBAdaptor *dba) {
  RepeatConsensusAdaptor *rca;

  if ((rca = (RepeatConsensusAdaptor *)calloc(1,sizeof(RepeatConsensusAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for RepeatConsensusAdaptor\n");
    return NULL;
  }
  BaseAdaptor_init((BaseAdaptor *)rca, dba, REPEATCONSENSUS_ADAPTOR);

  return rca;
}

RepeatConsensus *RepeatConsensusAdaptor_fetchByDbID(RepeatConsensusAdaptor *rca, int64 dbID) {

/*
    my ($rc) = @{$self->_generic_fetch("repeat_consensus_id = $db_id")}; 

    return $rc;   
*/
}

RepeatConsensus *RepeatConsensusAdaptor_fetchByName(RepeatConsensusAdaptor *rca, char *name) {

/*
    my ($rc) = @{$self->_generic_fetch("repeat_name = '$name'")};   

    return $rc;
*/
}

RepeatConsensus *RepeatConsensusAdaptor_fetchByNameAndClass(RepeatConsensusAdaptor *rca, char *name, char *class) {

/*
    my ($rc) = @{$self->_generic_fetch(qq{
      repeat_name  = '$name'
      AND repeat_class = '$class'
    })};

    return $rc;
*/
}


Set *RepeatConsensusAdaptor_fetchByClassAndSeq(RepeatConsensusAdaptor *rca, char *class, char *seq) {

/*
    return $self->_generic_fetch(qq{
            repeat_class     = '$class'
        AND repeat_consensus = '$seq'
    });
*/
}

Set *RepeatConsensusAdaptor_genericFetch(RepeatConsensusAdaptor *rca, char *whereClause) {
  StatementHandle *sth;
  ResultRow *row;
  char qStr[1024];
  Set *consensi;
    
  sprintf(qStr,"SELECT repeat_consensus_id, repeat_name,"
               "       repeat_class, LENGTH(repeat_consensus)"
               " FROM repeat_consensus"
               " WHERE %s", whereClause);

  rca->prepare((BaseAdaptor *)rca,qStr,strlen(qStr));
  sth->execute(sth);

  consensi = Set_new();
  
  while (row = sth->fetchRow(sth)) {
    RepeatConsensus *rc = RepeatConsensus_new();

    RepeatConsensus_setDbID(rc, row->getLongLongAt(row,0));
    RepeatConsensus_setName(rc, row->getStringAt(row,1));
    RepeatConsensus_setRepeatClass(rc, row->getStringAt(row,2));
    // Do I really need to do this???? RepeatConsensus_setLength(rc, row->getIntAt(row,3));
    RepeatConsensus_setAdaptor(rc, rca);

    Set_addElement(consensi, rc);
  }
  return consensi;
}

int RepeatConsensusAdaptor_store(RepeatConsensusAdaptor *rca, Set *consensi) {
/*
  my( $self, @consensi ) = @_;
  
  my $sth = $self->prepare(q{
    INSERT into repeat_consensus( repeat_consensus_id
          , repeat_name
          , repeat_class
          , repeat_consensus )
      VALUES (NULL, ?,?,?)
    });
    
  foreach my $rc (@consensi) {
    my $name  = $rc->name
      or $self->throw("name not set");
    my $class = $rc->repeat_class
      or $self->throw("repeat_class not set");
    my $seq   = $rc->repeat_consensus
      or $self->throw("repeat_consensus not set");
    
    $sth->execute($name, $class, $seq);
    
    my $db_id = $sth->{'mysql_insertid'}
    or $self->throw("Didn't get an insertid from the INSERT statement");
    
    $rc->dbID($db_id);
    $rc->adaptor($self);
  }
*/
}

