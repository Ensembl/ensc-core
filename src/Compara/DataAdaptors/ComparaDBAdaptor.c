#include "ComparaDBAdaptor.h"

#include "DBAdaptor.h"
#include "Species.h"
#include "MetaContainer.h"
#include "StrUtil.h"
#include "FileUtil.h"

Vector *ComparaDBAdaptor_readConfFile(ComparaDBAdaptor *cdba, char *fileName);

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
    //read configuration file from disk
    ComparaDBAdaptor_readConfFile(cdba, confFile);
  }
  return cdba;
}

/* Format of conf file:-
# $compara_db->get_bd_adaptor("Homo sapiens","NCBI31")

[
['Homo sapiens', 'NCBI31', {'host' => "ecs2b.internal.sanger.ac.uk",
                            'user' => "ensro",
                            'dbname' => "ens_NCBI_31",
                            'module' => 'Bio::EnsEMBL::DBSQL::DBAdaptor'}],
]
*/
Vector *ComparaDBAdaptor_readConfFile(ComparaDBAdaptor *cdba, char *fileName) {
  int inOuter = 0;
  int inInner = 0;
  int inHash = 0;
  FILE *confFP;
  char line[MAXSTRLEN];
  char *confStr;
  int confStrLen = 0;
  char *chP;

  if ((confFP = FileUtil_open(fileName, "r", "ComparaDBAdaptor_readConfFile")) == NULL) {
    fprintf(stderr, "Error: Failed opening conf file.\n");
    exit(1);
  }

  StrUtil_copyString(&confStr, "", 0);
  while (fgets(line, MAXSTRLEN, confFP)) {
    int len;

    StrUtil_truncateAtChar(line, '#');
    len = StrUtil_rmspace(line);

    if (line[0] == '[' || inOuter) {
      inOuter = 1;
      confStr = StrUtil_appendString(confStr, line);
      confStrLen+=len;
    } else if (len > 0) {
      fprintf(stderr,"Error: Failed reading conf file at %s\n",line);
    }
  }
  
  // printf("confStr = %s\n",confStr);

  if (confStr[confStrLen-1] != ']') {
    fprintf(stderr,"Error: conf file missing closing ']'\n");
    exit(1);
  }

// Remove the outer []
  chP = confStr+1;
  confStr[confStrLen-1] = '\0';

  while (*chP != '\0') {
    int sectLen;
    int tokLen;
    char *startSect;
    int tokNum = 0;
    int port = 3306;
    char species[1024];
    char assembly[1024];
    char user[1024];
    char dbname[1024];
    char host[1024];
    char genomeHashKey[1024];
    DBAdaptor *dba;

    host[0] = species[0] = assembly[0] = user[0] = dbname[0] = '\0';

    if (*chP != '[') {
      fprintf(stderr, "Error: didn't get a '[' at start of conf section\n");
      exit(1);
    }

    startSect = ++chP;

    sectLen = StrUtil_truncateAtChar(chP, ']');
    // printf("section = %s\n",chP);

    while ((tokLen = StrUtil_truncateAtChar(chP,','))) {
      //printf("token = %s\n",chP);
      switch (tokNum) {
        case 0:
          strcpy(species, chP);
          StrUtil_rmQuotes(species);
          break;
        case 1:
          strcpy(assembly, chP);
          StrUtil_rmQuotes(assembly);
          break;
        default:
          {
            char key[1024];
            char value[1024];
            char *fillP = key;
            char *tokP = chP;
            char prevCh = '\0';
  
            while (*tokP != '\0') {
              if (prevCh == '=' && *tokP == '>') {
                *(fillP-1) = '\0';
                fillP = value;
              } else if (*tokP!='{' && *tokP!='}') {
                *fillP = *tokP;
                fillP++;
              }
              prevCh = *tokP;
              tokP++;
            }
            *fillP = '\0';

            StrUtil_rmQuotes(key);
            StrUtil_rmQuotes(value);
            if (!strcmp(key,"host")) {
              strcpy(host, value);
            } else if (!strcmp(key,"user")) {
              strcpy(user, value);
            } else if (!strcmp(key,"port")) {
              port = atoi(value);
            } else if (!strcmp(key,"dbname")) {
              strcpy(dbname, value);
            } else if (!strcmp(key,"module")) {
              if (strcmp(value,"Bio::EnsEMBL::DBSQL::DBAdaptor")) {
                fprintf(stderr,"Error: Module not normal adaptor is %s\n",value);
                exit(1);
              }
            }
          }
          break;
      }

      tokNum++;  
      chP+=tokLen+1;
    }

    if (host[0] == '\0' ||
        species[0] == '\0' ||
        assembly[0] == '\0' ||
        user[0] == '\0' ||
        dbname[0] == '\0') {
      fprintf(stderr,"Error: Missing parameter in conf section\n");
      exit(1);
    }
    dba = DBAdaptor_new(host,user,NULL,dbname,port,NULL);
    DBAdaptor_setAssemblyType(dba, assembly);
    sprintf(genomeHashKey,"%s:%s",species,assembly);
    StringHash_add(cdba->genomes, genomeHashKey, dba);

    chP = startSect + (sectLen+2);
  }


  free(confStr);
  FileUtil_close(confFP,fileName);
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

/* Seems to be crap
SyntenyRegionAdaptor *ComparaDBAdaptor_getSyntenyRegionAdaptor(ComparaDBAdaptor *cdba) {
  if (!DBConnection_getAdaptor(cdba->dbc,SYNTENYREGION_ADAPTOR)) {
    DBConnection_addAdaptor(cdba->dbc,
                            (BaseAdaptor *)SyntenyRegionAdaptor_new(cdba));
  }
  return (SyntenyRegionAdaptor *)DBConnection_getAdaptor(cdba->dbc,SYNTENYREGION_ADAPTOR);
}
*/

ComparaDNAAlignFeatureAdaptor *ComparaDBAdaptor_getComparaDNAAlignFeatureAdaptor(ComparaDBAdaptor *cdba)  {
  if (!DBConnection_getAdaptor(cdba->dbc,COMPARADNAALIGNFEATURE_ADAPTOR)) {
    DBConnection_addAdaptor(cdba->dbc,
                            (BaseAdaptor *)ComparaDNAAlignFeatureAdaptor_new(cdba));
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
