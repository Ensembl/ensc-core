/*
 * Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *      http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "GenomeDBAdaptor.h"


GenomeDBAdaptor *GenomeDBAdaptor_new(ComparaDBAdaptor *dba) {
  GenomeDBAdaptor *gda;

  if ((gda = (GenomeDBAdaptor *)calloc(1,sizeof(GenomeDBAdaptor))) == NULL) {
    fprintf(stderr,"Error: Failed allocating gda\n");
    exit(1);
  }
  BaseComparaAdaptor_init((BaseComparaAdaptor *)gda, dba, GENOMEDB_ADAPTOR);


  return gda;
}

GenomeDB *GenomeDBAdaptor_fetchByDbID(GenomeDBAdaptor *gda, IDType dbID) {
  GenomeDB *gdb;
  DBAdaptor *dba;

  if (!dbID) {
    fprintf(stderr,"Error: Must have dbID to fetch by dbid\n");
    exit(1);
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
  dba = ComparaDBAdaptor_getDBAdaptor(gda->dba, GenomeDB_getName(gdb), 
                                      GenomeDB_getAssembly(gdb));

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

Vector *GenomeDBAdaptor_fetchAll(GenomeDBAdaptor *gda) {
  Vector *genomeDBs;
  IDType *keys;
  int i;
  
  if (!gda->genomeDBCache) {
    GenomeDBAdaptor_createGenomeDBs(gda);
  }

  keys = IDHash_getKeys(gda->genomeDBCache);
  genomeDBs = Vector_new();

  for (i=0; i<IDHash_getNumValues(gda->genomeDBCache); i++) {
    GenomeDB *gdb = IDHash_getValue(gda->genomeDBCache, keys[i]);
    ComparaDBAdaptor *cdb = gda->dba;
    DBAdaptor *dba = ComparaDBAdaptor_getDBAdaptor(cdb, GenomeDB_getName(gdb), 
                                                   GenomeDB_getAssembly(gdb));
    Vector_addElement(genomeDBs, gdb);
    if (dba) {
      GenomeDB_setDBAdaptor(gdb, dba);
    }
  }

  free(keys);
    
  return genomeDBs;
} 

GenomeDB *GenomeDBAdaptor_fetchByNameAssembly(GenomeDBAdaptor *gda, char *name, char *assembly) {
  StatementHandle *sth;
  ResultRow *row;
  char qStr[512];
  IDType id;

  if (!name || !assembly) {
    fprintf(stderr,"Error: name and assembly arguments are required\n");
    exit(1);
  }

  sprintf(qStr,
	     "SELECT genome_db_id"
             " FROM genome_db"
             " WHERE name = '%s' AND assembly = '%s'", name, assembly);
  sth = gda->prepare((BaseAdaptor *)gda, qStr, strlen(qStr));

  sth->execute(sth);

  if (!(row = sth->fetchRow(sth))) {
    fprintf(stderr,"Error: No GenomeDB with this name [%s] and " 
                   "assembly [%s]\n", name, assembly);
    exit(1);
  }

  id = row->getLongLongAt(row,0);

  return GenomeDBAdaptor_fetchByDbID(gda, id);
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
        GenomeDB_getName(gdb),GenomeDB_getAssembly(gdb), GenomeDB_getTaxonId(gdb));

    sth = gda->prepare((BaseAdaptor *)gda, qStr, strlen(qStr));
    sth->execute(sth);
    dbID = sth->getInsertId(sth);
    sth->finish(sth);
  }

  // update the genomeDB object so that it's dbID and adaptor are set
  GenomeDB_setDbID(gdb, dbID);
  GenomeDB_setAdaptor(gdb, (BaseAdaptor *)gda);

  return dbID;
}

void GenomeDBAdaptor_createGenomeDBs(GenomeDBAdaptor *gda) {
  char qStr[1024];
  StatementHandle *sth;
  ResultRow *row;

  // Populate the hash array which cross-references the consensus
  // and query dbs
 
  gda->genomeDBCache = IDHash_new(IDHASH_SMALL);
  gda->genomeConsensusXrefList = StringHash_new(STRINGHASH_SMALL);
  gda->genomeQueryXrefList = StringHash_new(STRINGHASH_SMALL);

  sprintf(qStr,
     "SELECT consensus_genome_db_id, query_genome_db_id, method_link_id"
     " FROM genomic_align_genome");

  sth = gda->prepare((BaseAdaptor *)gda, qStr, strlen(qStr));
  sth->execute(sth);

  while ((row = sth->fetchRow(sth))) {
    char cKey[512];
    char qKey[512];
    Vector *qVector;
    Vector *cVector;
    IDType *idP;
    IDType cons  = row->getLongLongAt(row,0);
    IDType query = row->getLongLongAt(row,1);
    IDType methId = row->getLongLongAt(row,2);

    sprintf(cKey,IDFMTSTR ":" IDFMTSTR, cons, methId);
    sprintf(qKey,IDFMTSTR ":" IDFMTSTR, query, methId);

    if (!StringHash_contains(gda->genomeConsensusXrefList, cKey)) {
      StringHash_add(gda->genomeConsensusXrefList, cKey, Vector_new());
    }
    if (!StringHash_contains(gda->genomeQueryXrefList, qKey)) {
      StringHash_add(gda->genomeQueryXrefList, qKey, Vector_new());
    }

    cVector =  StringHash_getValue(gda->genomeConsensusXrefList, cKey);
    if ((idP = (IDType *)calloc(1,sizeof(IDType))) == NULL) {
      fprintf(stderr,"Error: Failed allocating idP\n");
      exit(1);
    }
    *idP = query;
    Vector_addElement(cVector, idP);

    qVector =  StringHash_getValue(gda->genomeQueryXrefList, qKey);
    if ((idP = (IDType *)calloc(1,sizeof(IDType))) == NULL) {
      fprintf(stderr,"Error: Failed allocating idP\n");
      exit(1);
    }
    *idP = cons;
    Vector_addElement(qVector, idP);
  }

  sth->finish(sth);

  // grab all the possible species databases in the genome db table
  sprintf(qStr,
     "SELECT genome_db_id, name, assembly, taxon_id"
     " FROM genome_db");

  sth = gda->prepare((BaseAdaptor *)gda, qStr, strlen(qStr));
  sth->execute(sth);

  // build a genome db for each species
  while ((row = sth->fetchRow(sth))) {

    GenomeDB *gdb = GenomeDB_new();
    GenomeDB_setDbID(gdb, row->getLongLongAt(row,0));
    GenomeDB_setName(gdb, row->getStringAt(row,1));
    GenomeDB_setAssembly(gdb, row->getStringAt(row,2));
    GenomeDB_setTaxonId(gdb, row->getLongLongAt(row,3));
    GenomeDB_setAdaptor(gdb,(BaseAdaptor *)gda);

    IDHash_add(gda->genomeDBCache,GenomeDB_getDbID(gdb),gdb);
  }
  sth->finish(sth);
}

int GenomeDBAdaptor_checkForConsensusDb(GenomeDBAdaptor *gda, GenomeDB *queryGdb, 
                                        GenomeDB *conGdb, IDType methodLinkId) {
  char key[512];
  IDType cid = GenomeDB_getDbID(conGdb);
  IDType qid = GenomeDB_getDbID(queryGdb);
  
  sprintf(key,IDFMTSTR ":" IDFMTSTR, cid, methodLinkId);
  
  if (StringHash_contains(gda->genomeConsensusXrefList,key)) {
    int i;
    Vector *cxl = StringHash_getValue(gda->genomeConsensusXrefList, key);

    for (i=0; i<Vector_getNumElement(cxl); i++) {
      IDType *id = Vector_getElementAt(cxl, i);
      if (qid == *id) {
        return 1;
      }
    }
  }
  return 0;
}

int GenomeDBAdaptor_checkForQueryDb(GenomeDBAdaptor *gda, GenomeDB *conGdb, 
                                    GenomeDB *queryGdb, IDType methodLinkId) {

  char key[512];
  IDType cid = GenomeDB_getDbID(conGdb);
  IDType qid = GenomeDB_getDbID(queryGdb);
  
  sprintf(key,IDFMTSTR ":" IDFMTSTR, qid, methodLinkId);
  
  if (StringHash_contains(gda->genomeQueryXrefList,key)) {
    int i;
    Vector *qxl = StringHash_getValue(gda->genomeQueryXrefList, key);

    for (i=0; i<Vector_getNumElement(qxl); i++) {
      IDType *id = Vector_getElementAt(qxl, i);
      if (cid == *id) {
	return 1;
      }
    }
  }
  return 0;
}



/*
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
*/

Vector *GenomeDBAdaptor_getAllDbLinks(GenomeDBAdaptor *gda, GenomeDB *refGdb, IDType methodLinkId) {
  IDType id = GenomeDB_getDbID(refGdb);
  Vector *gdbList = Vector_new();
  char key[512];

  sprintf(key, IDFMTSTR ":" IDFMTSTR, id, methodLinkId);

  // check for occurences of the db we are interested in
  // in the consensus list of dbs

  if (StringHash_contains(gda->genomeConsensusXrefList, key)) {
    int i;
    Vector *cxl = StringHash_getValue(gda->genomeConsensusXrefList, key);

    for (i=0; i<Vector_getNumElement(cxl); i++) {
      IDType *id = Vector_getElementAt(cxl, i);
      Vector_addElement(gdbList, IDHash_getValue(gda->genomeDBCache, *id));
    }
  }

  // and check for occurences of the db we are interested in
  // in the query list of dbs
  if (StringHash_contains(gda->genomeQueryXrefList, key)) {
    int i;
    Vector *qxl = StringHash_getValue(gda->genomeQueryXrefList, key);

    for (i=0; i<Vector_getNumElement(qxl); i++) {
      IDType *id = Vector_getElementAt(qxl, i);
      Vector_addElement(gdbList, IDHash_getValue(gda->genomeDBCache, *id));
    }
  }

  return gdbList;
}

