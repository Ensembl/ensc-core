#include "GenomeDBAdaptor.h"


package Bio::EnsEMBL::Compara::DBSQL::GenomeDBAdaptor;
use vars qw(@ISA);
use strict;


use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Compara::GenomeDB;


# Hashes for storing a cross-referencing of compared genomes
my %genome_consensus_xreflist;
my %genome_query_xreflist;


@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

GenomeDBAdaptor *GenomeDBAdaptor_new() {
}

GenomeDB *GenomeDBAdaptor_fetchByDbID(GenomeDBAdaptor *gda, IDType dbID) {
  GenomeDB *gdb;
  DBAdaptor *dba;

  if (!dbID) {
    $self->throw("Must fetch by dbid");
  }

  // check to see whether all the GenomeDBs have already been created
  if (!gda->genomeDBCache) {
    GenomeDBAdaptor_createGenomeDBs(gda);
  }

  if (!IDHash_contains(gda->genomeDBCache, dbID)) {
    return NULL; // return undef if fed a bogus dbID
  }
  gdb = IDHash_getValue(gda->genomeDBCache, dbID);

  // set up the dbadaptor for this genome db
  // this could have been added after the cache was created which is why
  // it is re-added every request
  my $dba = $self->db->get_db_adaptor($gdb->name, $gdb->assembly);

  if (!dba) {
    fprintf(stderr,"Error: Could not obtain DBAdaptor for dbID [" IDFMTSTR "].\n" 
                   "Genome DBAdaptor for name=[%s], "
                   "assembly=[%s] must be loaded using config file or\n" 
                   "ComparaDBAdaptor_add_genome\n", 
                  dbID, GenomeDB_getName(gdb), GenomeDB_getAssembly(gdb));
    exit(1);
  }

  GenomeDB_setDBAdaptor(gdb, dba);

  return gdb;
}


=head2 fetch_all

  Args       : none
  Example    : none
  Description: gets all GenomeDBs for this compara database
  Returntype : listref Bio::EnsEMBL::Compara::GenomeDB
  Exceptions : none
  Caller     : general

=cut

Vector *GenomeDBAdaptor_fetchAll(GenomeDBAdaptor *gda) {

  if ( !defined $self->{'_GenomeDB_cache'}) {
    $self->create_GenomeDBs;
  }

  my @genomeDBs = values %{$self->{'_cache'}};

  for my $gdb ( @genomeDBs ) {
    my $dba = $self->db->get_db_adaptor($gdb->name, $gdb->assembly);
    if($dba) {
      $gdb->db_adaptor($dba);
    }
  }
    
  return genomeDBs;
} 



=head2 fetch_by_name_assembly

  Arg [1]    : string $name
  Arg [2]    : string $assembly
  Example    : $gdb = $gdba->fetch_by_name_assembly("Homo sapiens", 'NCBI_31');
  Description: Retrieves a genome db using the name of the species and
               the assembly version.
  Returntype : Bio::EnsEMBL::Compara::GenomeDB
  Exceptions : thrown if GenomeDB of name $name and $assembly cannot be found
  Caller     : general

=cut

GenomeDB *GenomeDBAdaptor_fetchByNameAssembly(GenomeDBAdaptor *gda, char *name, char *assembly) {

   unless($name && $assembly) {
     $self->throw('name and assembly arguments are required');
   }

   my $sth = $self->prepare(
	     "SELECT genome_db_id
              FROM genome_db
              WHERE name = ? AND assembly = ?");

   $sth->execute($name, $assembly);

   my ($id) = $sth->fetchrow_array();

   if( !defined $id ) {
       $self->throw("No GenomeDB with this name [$name] and " .
		    "assembly [$assembly]");
   }

   return $self->fetch_by_dbID($id);
}

IDType GenomeDBAdaptor_store(GenomeDBAdaptor *gda, GenomeDB *gdb) {
  IDType dbID = 0;
  StatementHandle *sth;
  ResultRow *row;
  char qStr[1024];

  if (!gdb) {
    fprintf(stderr, "Error: Must have genomedb arg in store\n");
    exit(1);
  }

  if (!GenomeDB_getName(gdb) || !GenomeDB_getAssembly(gdb) || !GenomeDB_getTaxonId(gdb)) {
    fprintf(stderr, "Error: genome db must have a name, assembly, and taxon_id\n");
    exit(1);
  }

  sprintf(qStr, 
      "SELECT genome_db_id"
      " FROM genome_db"
      " WHERE name = '%s' and assembly = '%s'", 
      GenomeDB_getName(gdb), GenomeDB_getAssembly(gdb));

  sth = gda->prepare((BaseAdaptor *)gda, qStr,strlen(qStr));
  sth->execute(sth);

  if ((row = sth->fetchRow(sth))) {
    dbID = row->getLongLongAt(row,0);
  }
  sth->finish(sth);

  if (!dbID) {
    // if the genome db has not been stored before, store it now
    sprintf(qStr,
        "INSERT into genome_db (name,assembly,taxon_id)"
        " VALUES ('%s','%s', %d)",
        GenomeDB_getName(gdb),GenomeDB_getAssemby(gdb), GenomeDB_getTaxonId(gdb));

    sth = gda->prepare((BaseAdaptor *)gda, qStr, strlen(qStr));
    sth->execute(sth);
    dbID = sth->getInsertId(sth);
  }

  // update the genomeDB object so that it's dbID and adaptor are set
  GenomeDB_setDbID(gdb, dbID);
  GenomeDB_setAdaptor(gdb, dfa);

  return dbID;
}

void GenomeDBAdaptor_createGenomeDBs(GenomeDBAdaptor *gda) {

  // Populate the hash array which cross-references the consensus
  // and query dbs

  my $sth = $self->prepare("
     SELECT consensus_genome_db_id, query_genome_db_id, method_link_id
     FROM genomic_align_genome
  ");

  $sth->execute;

  while ( my @db_row = $sth->fetchrow_array() ) {
    my ( $con, $query, $method_link_id ) = @db_row;

    $genome_consensus_xreflist{$con .":" .$method_link_id} ||= [];
    $genome_query_xreflist{$query .":" .$method_link_id} ||= [];

    push @{ $genome_consensus_xreflist{$con .":" .$method_link_id}}, $query;
    push @{ $genome_query_xreflist{$query .":" .$method_link_id}}, $con;
  }

  // grab all the possible species databases in the genome db table
  $sth = $self->prepare("
     SELECT genome_db_id, name, assembly, taxon_id
     FROM genome_db 
   ");
   $sth->execute;

  // build a genome db for each species
  while ( my @db_row = $sth->fetchrow_array() ) {
    my ($dbid, $name, $assembly, $taxon_id) = @db_row;

    GenomeDB *gdb = GenomeDB_new();
    GenomeDB_setDbID(gdb, row->getLongLongAt(row,0));
    GenomeDB_setName(gdb, row->getStringAt(row,1));
    $gdb->assembly($assembly);
    $gdb->taxon_id($taxon_id);
    GenomeDB_setAdaptor(gdb,gda);

    IDHash_add(gda->genomeDBCache,GenomeDB_getDbID(gdb),gdb);
  }

}


=head2 check_for_consensus_db

  Arg[1]     : Bio::EnsEMBL::Compara::GenomeDB $consensus_genomedb
  Arg[2]     : Bio::EnsEMBL::Compara::GenomeDB $query_genomedb
  Arg[3]     : int $method_link_id
  Example    :
  Description: Checks to see whether a consensus genome database has been
               analysed against the specific query genome database.
               Returns the dbID of the database of the query genomeDB if 
               one is found.  A 0 is returned if no match is found.
  Returntype : int ( 0 or 1 )
  Exceptions : none
  Caller     : Bio::EnsEMBL::Compara::GenomeDB.pm

=cut


int GenomeDBAdaptor_checkForConsensusDb(GenomeDBAdaptor *gda, GenomeDB *queryGdb, GenomeDB *conGdb, IDType methodLinkId) {

  // just to make things a wee bit more readable
  my $cid = $con_gdb->dbID;
  my $qid = $query_gdb->dbID;
  
  if ( exists $genome_consensus_xreflist{$cid .":" .$method_link_id} ) {
    for my $i ( 0 .. $#{$genome_consensus_xreflist{$cid .":" .$method_link_id}} ) {
      if ( $qid == $genome_consensus_xreflist{$cid .":" .$method_link_id}[$i] ) {
	return 1;
      }
    }
  }
  return 0;
}


=head2 check_for_query_db

  Arg[1]     : Bio::EnsEMBL::Compara::GenomeDB $query_genomedb
  Arg[2]     : Bio::EnsEMBL::Compara::GenomeDB $consensus_genomedb
  Arg[3]     : int $method_link_id
  Example    : none
  Description: Checks to see whether a query genome database has been
               analysed against the specific consensus genome database.
               Returns the dbID of the database of the consensus 
               genomeDB if one is found.  A 0 is returned if no match is
               found.
  Returntype : int ( 0 or 1 )
  Exceptions : none
  Caller     : Bio::EnsEMBL::Compara::GenomeDB.pm

=cut

int GenomeDBAdaptor_checkForQueryDb(GenomeDBAdaptor *gda, GenomeDB *conGdb, GenomeDB *queryGdb, IDType methodLinkId) {

  // just to make things a wee bit more readable
  my $cid = $con_gdb->dbID;
  my $qid = $query_gdb->dbID;

  if ( exists $genome_query_xreflist{$qid .":" .$method_link_id} ) {
    for my $i ( 0 .. $#{$genome_query_xreflist{$qid .":" .$method_link_id}} ) {
      if ( $cid == $genome_query_xreflist{$qid .":" .$method_link_id}[$i] ) {
	return 1;
      }
    }
  }
  return 0;
}



=head2 get_all_db_links

  Arg[1]     : Bio::EnsEMBL::Compara::GenomeDB $query_genomedb
  Arg[2]     : int $method_link_id
  Example    : 
  Description: For the GenomeDB object passed in, check is run to
               verify which other genomes it has been analysed against
               irrespective as to whether this was as the consensus
               or query genome. Returns a list of matching dbIDs 
               separated by white spaces. 
  Returntype : listref of Bio::EnsEMBL::Compara::GenomeDBs 
  Exceptions : none
  Caller     : Bio::EnsEMBL::Compara::GenomeDB.pm

=cut

sub GenomeDBAdaptor_getAllDbLinks(GenomeDBAdaptor *gda, GenomeDB *refGdb, IDType methodLinkId) {
  my ( $self, $ref_gdb,$method_link_id ) = @_;
  
  my $id = $ref_gdb->dbID;
  my @gdb_list;

  // check for occurences of the db we are interested in
  // in the consensus list of dbs
  if ( exists $genome_consensus_xreflist{$id . ":" .$method_link_id} ) {
    for my $i ( 0 .. $#{ $genome_consensus_xreflist{$id . ":" .$method_link_id} } ) {
      push @gdb_list, $self->{'_cache'}->{$genome_consensus_xreflist{$id . ":" .$method_link_id}[$i]};
    }
  }

  // and check for occurences of the db we are interested in
  // in the query list of dbs
  if ( exists $genome_query_xreflist{$id . ":" .$method_link_id} ) {
    for my $i ( 0 .. $#{ $genome_query_xreflist{$id . ":" .$method_link_id} } ) {
      push @gdb_list, $self->{'_cache'}->{$genome_query_xreflist{$id . ":" .$method_link_id}[$i]};
    }
  }

  return \@gdb_list;
}

