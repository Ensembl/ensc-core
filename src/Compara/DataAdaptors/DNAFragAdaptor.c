/*
 * Copyright [1999-2017] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

#include "DNAFragAdaptor.h"
#include "DNAFrag.h"
#include "Vector.h"
#include "GenomeDBAdaptor.h"

DNAFragAdaptor *DNAFragAdaptor_new(ComparaDBAdaptor *dba) {
  DNAFragAdaptor *dfa = NULL;

  if ((dfa = (DNAFragAdaptor *)calloc(1,sizeof(DNAFragAdaptor))) == NULL) {
    fprintf(stderr,"Error: Failed allocating dfa\n");
    dfa = NULL;
  } else {
    BaseComparaAdaptor_init((BaseComparaAdaptor *)dfa, dba, DNAFRAG_ADAPTOR);
  }

  return dfa;
}

DNAFrag *DNAFragAdaptor_fetchByDbID(DNAFragAdaptor *dfa, IDType dbID) {
  char qStr[1024];
  Vector *dnaFrags;
  DNAFrag *firstFrag = NULL;
  StatementHandle *sth;

  if (!dbID) {
    fprintf(stderr,"Error: Must have a dbid to fetch by dbid\n");
  } else {
    sprintf(qStr,
            "SELECT genome_db_id, dnafrag_type, dnafrag_id,"
            "       name, start, end"
            "  FROM dnafrag"
            " WHERE dnafrag_id = " IDFMTSTR " limit 1", dbID);

    sth = dfa->prepare((BaseAdaptor *)dfa,qStr,strlen(qStr));
    sth->execute(sth);

    dnaFrags = DNAFragAdaptor_objectsFromStatementHandle(dfa,sth);

    if (!Vector_getNumElement(dnaFrags)) {
      fprintf(stderr,"Error: No dnafrag with dbID " IDFMTSTR "\n",dbID);
    } else {
      firstFrag = Vector_getElementAt(dnaFrags,0);
    }

    Vector_free(dnaFrags);
  }
  return firstFrag;
}

Vector *DNAFragAdaptor_fetchAllByGenomeDBRegion(DNAFragAdaptor *dfa, GenomeDB *genomeDB,
                      char *dnaFragType, char *name, int *startP, int *endP) {
  Vector *result = NULL;
  int ok = 1;
  IDType gdbId;
  char qStr[1024];
  StatementHandle *sth;

  if (!genomeDB) {
    fprintf(stderr, "Error: genome_db arg must be set\n");
    ok = 0;
  }

  if (ok) {
    gdbId = GenomeDB_getDbID(genomeDB);

    if (!gdbId) {
      fprintf(stderr, "Error: GenomeDB does not have a dbID. Is it stored in the db?\n");
      ok = 0;
    }
  }

  if (ok && !dnaFragType) {
    fprintf(stderr,"Error: dnafrag_type argument must be defined\n");
    ok = 0;
  }

  if (ok) {
    sprintf(qStr,
            "SELECT d.genome_db_id, d.dnafrag_type, d.dnafrag_id,"
            "        d.name, d.start, d.end"
            " FROM  dnafrag d"
            " WHERE d.dnafrag_type = '%s'"
            " AND   d.genome_db_id = " IDFMTSTR, dnaFragType, gdbId);


    if (name) {
      char tmpStr[128];
      sprintf(tmpStr," AND d.name = '%s'", name);
      strcat(qStr, tmpStr);
    }

    if (startP) {
      char tmpStr[64];
      sprintf(tmpStr," AND d.end >= %d", *startP);
      strcat(qStr, tmpStr);
    }

    if (endP) {
      char tmpStr[64];
      sprintf(tmpStr," AND d.start <= %d", *endP);
      strcat(qStr, tmpStr);
    }

    sth = dfa->prepare((BaseAdaptor *)dfa,qStr,strlen(qStr));
    sth->execute(sth);

    result = DNAFragAdaptor_objectsFromStatementHandle(dfa, sth);
  }

  return result;
}

Vector *DNAFragAdaptor_fetchAll(DNAFragAdaptor *dfa) {
  char qStr[256];
  StatementHandle *sth;
 
  sprintf(qStr,
           "SELECT genome_db_id, dnafrag_type, dnafrag_id," 
           "       name, start, end"
           " FROM dnafrag");
  sth = dfa->prepare((BaseAdaptor *)dfa,qStr,strlen(qStr));
  sth->execute(sth);

  return DNAFragAdaptor_objectsFromStatementHandle(dfa, sth);
}


Vector *DNAFragAdaptor_objectsFromStatementHandle(DNAFragAdaptor *dfa, StatementHandle *sth) {
  Vector *results = Vector_new();
  ResultRow *row;

  GenomeDBAdaptor *gda = ComparaDBAdaptor_getGenomeDBAdaptor(dfa->dba);

  while ((row = sth->fetchRow(sth))) {
    DNAFrag *dnaFrag = DNAFrag_new();

    DNAFrag_setDbID(dnaFrag, row->getLongLongAt(row,2));
    DNAFrag_setName(dnaFrag,row->getStringAt(row,3));
    DNAFrag_setType(dnaFrag,row->getStringAt(row,1));
    DNAFrag_setStart(dnaFrag,row->getIntAt(row,4));
    DNAFrag_setEnd(dnaFrag,row->getIntAt(row,5));
    DNAFrag_setGenomeDB(dnaFrag,GenomeDBAdaptor_fetchByDbID(gda, row->getLongLongAt(row,0)));

    Vector_addElement(results, dnaFrag);
  }

  sth->finish(sth);

  return results;
}

IDType DNAFragAdaptor_store(DNAFragAdaptor *dfa, DNAFrag *dnaFrag) {
  IDType result = 0;
  int ok = 1;
  GenomeDB *gdb;
  IDType gid;
  char qStr[512];
  StatementHandle *sth;

  if (!dnaFrag) {
    fprintf(stderr,"Error: Must have a DNAFrag object to store it\n");
    ok = 0;
  }

  if (ok && DNAFrag_getAdaptor(dnaFrag) == (BaseAdaptor *)dfa) {
    result = DNAFrag_getDbID(dnaFrag);
  } else if (ok) {
    gdb = DNAFrag_getGenomeDB(dnaFrag);

    if (!gdb) {
      fprintf(stderr, "Error: Must have genomedb attached to the dnafrag to store the dnafrag\n");
      ok = 0;
    }

    if (ok && !GenomeDB_getDbID(gdb)) {
      fprintf(stderr, "Error: Genomedb must be stored (no dbID). Store genomedb first\n");
      ok = 0;
    }

    if (ok && !DNAFrag_getName(dnaFrag)) {
      fprintf(stderr, "Error: dna frag must have a name\n");
      ok = 0;
    }
   
    if (ok) {
      gid = GenomeDB_getDbID(gdb);

      sprintf(qStr,
              "INSERT INTO dnafrag ( genome_db_id, dnafrag_type,"
              "                      name, start, end )"
              " VALUES (" IDFMTSTR ",'%s','%s',%d,%d)",
              gid, 
              DNAFrag_getType(dnaFrag) ? DNAFrag_getType(dnaFrag) : "\\N",
              DNAFrag_getName(dnaFrag),
              DNAFrag_getStart(dnaFrag),
              DNAFrag_getEnd(dnaFrag));

      sth = dfa->prepare((BaseAdaptor *)dfa,qStr,strlen(qStr));
      sth->execute(sth);

      DNAFrag_setDbID(dnaFrag, sth->getInsertId(sth));
      DNAFrag_setAdaptor(dnaFrag, (BaseAdaptor *)dfa);

      result = DNAFrag_getDbID(dnaFrag);
    }
  }

  return result;
}

IDType DNAFragAdaptor_storeIfNeeded(DNAFragAdaptor *dfa, DNAFrag *dnaFrag) {
  IDType result = 0;
  int ok = 1;
  GenomeDB *gdb;
  IDType gid;
  IDType dnaFragId = 0;
  StatementHandle *sth;
  ResultRow *row;
  char qStr[1024];

  if (!dnaFrag) {
    fprintf(stderr,"Error: Must have a DNAFrag object to store it\n");
    ok = 0;
  }

  if (ok && DNAFrag_getAdaptor(dnaFrag) == (BaseAdaptor *)dfa) {
    result =DNAFrag_getDbID(dnaFrag);
  } else if (ok) {
    gdb = DNAFrag_getGenomeDB(dnaFrag);
  
    if (!gdb) {
      fprintf(stderr, "Error: Must have genomedb attached to the dnafrag to store the dnafrag\n");
      ok = 0;
    }


    if (ok && !GenomeDB_getDbID(gdb)) {
      fprintf(stderr, "Error: Genomedb must be stored (no dbID). Store genomedb first\n");
      ok = 0;
    }

    if (ok && !DNAFrag_getName(dnaFrag)) {
      fprintf(stderr, "Error: dna frag must have a name\n");
      ok = 0;
    }
   
    if (ok) {
      gid = GenomeDB_getDbID(gdb);
      sprintf(qStr,
              "SELECT dnafrag_id"
              "  FROM dnafrag"
              " WHERE name= '%s'"
              "   AND genome_db_id= " IDFMTSTR
              "   AND start = %d"
              "   AND end = %d", 
              DNAFrag_getName(dnaFrag), gid, DNAFrag_getStart(dnaFrag), DNAFrag_getEnd(dnaFrag));

      sth = dfa->prepare((BaseAdaptor *)dfa, qStr, strlen(qStr));
      sth->execute(sth);

      if ((row = sth->fetchRow(sth))) {
        dnaFragId = row->getLongLongAt(row,0);
      }
      sth->finish(sth);

      if (dnaFragId) {
        // dnafrag already stored
        DNAFrag_setDbID(dnaFrag, dnaFragId);
        DNAFrag_setAdaptor(dnaFrag, (BaseAdaptor *)dfa);
        result = dnaFragId;
      } else {
        result =DNAFragAdaptor_store(dfa, dnaFrag);
      }
    }
  }

  return result;
}
