#include "ComparaDBAdaptor.h"

=head2 new

  Arg [..]   : list of named arguments.  See Bio::EnsEMBL::DBConnection.
               [-CONF_FILE] optional name of a file containing configuration
               information for comparas genome databases.  If databases are
               not added in this way, then they should be added via the
               method add_DBAdaptor. An example of the conf file can be found
               in ensembl-compara/modules/Bio/EnsEMBL/Compara/Compara.conf.example
  Example    :  $db = new Bio::EnsEMBL::Compara::DBSQL::DBAdaptor(
						    -user   => 'root',
						    -dbname => 'pog',
						    -host   => 'caldy',
						    -driver => 'mysql',
                                                    -conf_file => 'conf.pl');
  Description: Creates a new instance of a DBAdaptor for the compara database.
  Returntype : Bio::EnsEMBL::Compara::DBSQL::DBAdaptor
  Exceptions : none
  Caller     : general

=cut

ComparaDBAdaptor *ComparaDBAdaptor_new(char *host, char *user, char *pass, char *dbname,
                                       unsigned int port, char *confFile) {
  ComparaDBAdaptor *cdba;

  if ((cdba = (ComparaDBAdaptor *)calloc(1,sizeof(ComparaDBAdaptor))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating ComparaDBAdaptor\n");
    return NULL;
  }

  cdba->dbc = DBConnection_new(host,user,pass,dbname,port);

  $self->{'genomes'} = {};
  if (confFile) {
    fprintf(stderr,"Error: Conf file reading not implemented\n");
    exit(1);

    //read configuration file from disk
    my @conf = @{do $conf_file};

    foreach my $genome (@conf) {
      my ($species, $assembly, $db_hash) = @$genome;
      my $db;

      my $module = $db_hash->{'module'};
      my $mod = $module;

      eval {
	// require needs /'s rather than colons
	if ( $mod =~ /::/ ) {
	  $mod =~ s/::/\//g;
	}
	require "${mod}.pm";

	$db = $module->new(-dbname => $db_hash->{'dbname'},
			   -host   => $db_hash->{'host'},
			   -user   => $db_hash->{'user'},
			   -pass   => $db_hash->{'pass'},
			   -port   => $db_hash->{'port'},
			   -driver => $db_hash->{'driver'});
      };

      if($@) {
	$self->throw("could not load module specified in configuration " .
		     "file:$@");
      }

      unless($db && ref $db && $db->isa('Bio::EnsEMBL::DBSQL::DBConnection')) {
	$self->throw("[$db] specified in conf file is not a " .
		     "Bio::EnsEMBL::DBSQL::DBConnection");
      }

      $self->{'genomes'}->{"$species:$assembly"} = $db;
    }
  }
  return cdba;
}


=head2 add_db_adaptor

  Arg [1]    : Bio::EnsEMBL::DBSQL::DBConnection
  Example    : $compara_db->add_db_adaptor($homo_sapiens_db);
  Description: Adds a genome-containing database to compara.  This database
               can be used by compara to obtain sequence for a genome on
               on which comparative analysis has been performed.  The database
               adaptor argument must define the get_MetaContainer argument
               so that species name and assembly type information can be
               extracted from the database.
  Returntype : none
  Exceptions : Thrown if the argument is not a Bio::EnsEMBL::DBConnection
               or if the argument does not implement a get_MetaContainer
               method.
  Caller     : general

=cut

sub ComparaDBAdaptor_addDBAdaptor(ComparaDBAdaptor *cdba, DBAdaptor *dba) {
  char key[1024];
  char *species;
  char *assembly;
  MetaContainer *mc;

  if (!dba) {
    fprintf(stderr, "Error: dba argument must be non null\n");
    exit(1):
  }

  mc = DBAdaptor_getMetaContainer(dba);

  species = $mc->get_Species->binomial;
  assembly = MetaContainer_getDefaultAssembly(mc);

  sprintf(key,"%s:%s",species,assembly);
  StringHash_add(cdba->genomes, key, dba);
}



=head2 get_db_adaptor

  Arg [1]    : string $species
               the name of the species to obtain a genome DBAdaptor for.
  Arg [2]    : string $assembly
               the name of the assembly to obtain a genome DBAdaptor for.
  Example    : $hs_db = $db->get_db_adaptor('Homo sapiens','NCBI_30');
  Description: Obtains a DBAdaptor for the requested genome if it has been
               specified in the configuration file passed into this objects
               constructor, or subsequently added using the add_genome
               method.  If the DBAdaptor is not available (i.e. has not
               been specified by one of the abbove methods) undef is returned.
  Returntype : Bio::EnsEMBL::DBSQL::DBConnection
  Exceptions : none
  Caller     : Bio::EnsEMBL::Compara::GenomeDBAdaptor

=cut

sub ComparaDBAdaptor_getDBAdaptor(ComparaDBAdaptor *cdba, char *species, char *assembly) {
  char key[1024];

  if (!species || !assembly) {
    fprintf(stderr, "Error: species and assembly arguments are required\n");
    exit(1);
  }

  sprintf(key,"%s:%s",species,assembly);
  if (!StringHash_contains(cdba->genomes, key)) {
    fprintf(stderr, "Error: No DBAdaptor loaded for %s\n", key);
    exit(1);
  }
  return StringHash_getValue(cdba->genomes, key);
}

SyntenyAdaptor *ComparaDBAdaptor_getSyntenyAdaptor(ComparaDBAdaptor *cdba) {
  if (!DBConnection_getAdaptor(cdba->dbc,SYNTENY_ADAPTOR)) {
    DBConnection_addAdaptor(cdba->dbc,
                            (BaseAdaptor *)SyntenyAdaptor_new(cdba));
  }
  return (SyntenyAdaptor *)DBConnection_getAdaptor(cdba->dbc,SYNTENY_ADAPTOR);
}

GenomeDBAdaptor *ComparaDBAdaptor_getGenomeDBAdaptor(ComparaDBAdaptor *cdba) {
  if (!DBConnection_getAdaptor(cdba->dbc,GENOMEDB_ADAPTOR)) {
    DBConnection_addAdaptor(cdba->dbc,
                            (BaseAdaptor *)GenomeDBAdaptor_new(cdba));
  }
  return (GenomeDBAdaptor *)DBConnection_getAdaptor(cdba->dbc,GENOMEDB_ADAPTOR);
}

DNAFragAdaptor *ComparaDBAdaptor_getDNAFragAdaptor(ComparaDBAdaptor *cdba) {
  if (!DBConnection_getAdaptor(cdba->dbc,DNAFRAG_ADAPTOR)) {
    DBConnection_addAdaptor(cdba->dbc,
                            (BaseAdaptor *)DNAFragAdaptor_new(cdba));
  }
  return (DNAFragAdaptor *)DBConnection_getAdaptor(cdba->dbc,DNAFRAG_ADAPTOR);
}

GenomicAlignAdaptor *ComparaDBAdaptor_getGenomicAlignAdaptor(ComparaDBAdaptor *cdba) {
  if (!DBConnection_getAdaptor(cdba->dbc,GENOMICALIGN_ADAPTOR)) {
    DBConnection_addAdaptor(cdba->dbc,
                            (BaseAdaptor *)GenomicAlignAdaptor_new(cdba));
  }
  return (GenomicAlignAdaptor *)DBConnection_getAdaptor(cdba->dbc,GENOMICALIGN_ADAPTOR);
}

HomologyAdaptor *ComparaDBAdaptor_getHomologyAdaptor(ComparaDBAdaptor *cdba) {
  if (!DBConnection_getAdaptor(cdba->dbc,HOMOLOGY_ADAPTOR)) {
    DBConnection_addAdaptor(cdba->dbc,
                            (BaseAdaptor *)HomologyAdaptor_new(cdba));
  }
  return (HomologyAdaptor *)DBConnection_getAdaptor(cdba->dbc,HOMOLOGY_ADAPTOR);
}

SyntenyRegionAdaptor *ComparaDBAdaptor_getSyntenyRegionAdaptor(ComparaDBAdaptor *cdba) {
  if (!DBConnection_getAdaptor(cdba->dbc,SYNTENYREGION_ADAPTOR)) {
    DBConnection_addAdaptor(cdba->dbc,
                            (BaseAdaptor *)SyntenyRegionAdaptor_new(cdba));
  }
  return (SyntenyRegionAdaptor *)DBConnection_getAdaptor(cdba->dbc,SYNTENYREGION_ADAPTOR);
}


ComparaDNAAlignFeatureAdaptor *ComparaDBAdaptor_getComparaDNAAlignFeatureAdaptor(ComparaDBAdaptor *cdba)  {
  if (!DBConnection_getAdaptor(cdba->dbc,COMPARADNAALIGNFEATURE_ADAPTOR)) {
    DBConnection_addAdaptor(cdba->dbc,
                            (BaseAdaptor *)ComparaDNAFeatureAdaptor_new(cdba));
  }
  return (ComparaDNAFeatureAdaptor *)DBConnection_getAdaptor(cdba->dbc,COMPARADNAALIGNFEATURE_ADAPTOR);
}

MetaContainer *ComparaDBAdaptor_getMetaContainer(ComparaDBAdaptor *cdba)  {
  if (!cdba->metaContainer) {
    cdba->metaContainer = MetaContainer_new(cdba);
  }
  return cdba->metaContainer;
}
