#include "ComparaDBAdaptor.h"

#include "DBAdaptor.h"
#include "Species.h"
#include "MetaContainer.h"

ComparaDBAdaptor *ComparaDBAdaptor_new(char *host, char *user, char *pass, char *dbname,
                                       unsigned int port, char *confFile) {
  ComparaDBAdaptor *cdba;

  if ((cdba = (ComparaDBAdaptor *)calloc(1,sizeof(ComparaDBAdaptor))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating ComparaDBAdaptor\n");
    return NULL;
  }

  cdba->dbc = DBConnection_new(host,user,pass,dbname,port);

  cdba->genomes = StringHash_new(STRINGHASH_SMALL);

  if (confFile) {
    fprintf(stderr,"Error: Conf file reading not implemented\n");
    exit(1);

#ifdef DONE
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
#endif
  }
  return cdba;
}

void ComparaDBAdaptor_addDBAdaptor(ComparaDBAdaptor *cdba, DBAdaptor *dba) {
  char key[1024];
  char *assembly;
  MetaContainer *mc;
  Species *species;
  char *speciesName;

  if (!dba) {
    fprintf(stderr, "Error: dba argument must be non null\n");
    exit(1);
  }

  mc = DBAdaptor_getMetaContainer(dba);

  species = MetaContainer_getSpecies(mc);
  speciesName = Species_getBinomialName(species);
  assembly = MetaContainer_getDefaultAssembly(mc);


  sprintf(key,"%s:%s",speciesName,assembly);
  StringHash_add(cdba->genomes, key, dba);

  Species_free(species);
// NIY ??
  free(assembly);
}

// Perl returns undef if not found - I throw
DBAdaptor *ComparaDBAdaptor_getDBAdaptor(ComparaDBAdaptor *cdba, char *species, char *assembly) {
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
  return (ComparaDNAAlignFeatureAdaptor *)DBConnection_getAdaptor(cdba->dbc,COMPARADNAALIGNFEATURE_ADAPTOR);
}

MetaContainer *ComparaDBAdaptor_getMetaContainer(ComparaDBAdaptor *cdba)  {
// HACK HACK HACK cast to DBAdaptor is BAD
  if (!cdba->metaContainer) {
    cdba->metaContainer = MetaContainer_new((DBAdaptor *)cdba);
  }
  return cdba->metaContainer;
}
