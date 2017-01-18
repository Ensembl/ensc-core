/*
 * Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
 * Copyright [2016-2017] EMBL-European Bioinformatics Institute
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

#include "HomologyAdaptor.h"
#include "StrUtil.h"


HomologyAdaptor *HomologyAdaptor_new(ComparaDBAdaptor *dba) {
  HomologyAdaptor *ha = NULL;

  if ((ha = (HomologyAdaptor *)calloc(1,sizeof(HomologyAdaptor))) == NULL) {
    fprintf(stderr,"Error: Failed allocating ha\n");
  } else {
    BaseComparaAdaptor_init((BaseComparaAdaptor *)ha, dba, HOMOLOGY_ADAPTOR);
  }

  return ha;
}


Vector *HomologyAdaptor_fetchHomologuesOfGeneInSpecies(HomologyAdaptor *ha, 
                          char *sp, char *gene, char *hSp) {
  char qStr[1024];
  Vector *genes;
  IDType *relationshipIds;
  int nRelationship;
  int i;
  char *hSpecies;
  char *species;

  species = StrUtil_copyString(&species,sp,0);
  species = StrUtil_strReplChr(species,'_',' ');

  hSpecies = StrUtil_copyString(&hSpecies,hSp,0);
  hSpecies = StrUtil_strReplChr(hSpecies,'_',' ');

  sprintf(qStr,
           "select grm.gene_relationship_id "
           " from   gene_relationship_member grm, "
           "        genome_db gd "
           " where  gd.genome_db_id = grm.genome_db_id "
           " and    gd.name = '%s' "
           " and    grm.member_stable_id = '%s' "
           " group by grm.gene_relationship_id", species, gene);

  nRelationship = HomologyAdaptor_getRelationships(ha,qStr,&relationshipIds);

  genes = Vector_new();
  for (i=0;i<nRelationship;i++) {
    Vector *homols = HomologyAdaptor_fetchHomologuesBySpeciesRelationshipId(ha,hSpecies,relationshipIds[i]);
    Vector_append(genes, homols);
    Vector_free(homols);
  }

  free(relationshipIds);
  free(species);
  free(hSpecies);

  return genes;
}


Vector *HomologyAdaptor_fetchHomologuesOfGene(HomologyAdaptor *ha, char *sp, char *gene) {
  char qStr[1024];
  char *species;
  int nRelationship;
  IDType *relationshipIds;
  Vector *genes;
  int i;

  species = StrUtil_copyString(&species,sp,0);
  species = StrUtil_strReplChr(species,'_',' ');

  sprintf(qStr,
            "select grm.gene_relationship_id "
            " from   gene_relationship_member grm, "
            "        genome_db gd "
            " where  gd.genome_db_id = grm.genome_db_id "
            " and    gd.name = '%s' "
            " and    grm.member_stable_id = '%s' "
            " group by grm.gene_relationship_id", species, gene);

  nRelationship = HomologyAdaptor_getRelationships(ha, qStr,&relationshipIds);

  genes = Vector_new();
  for (i=0;i<nRelationship;i++) {
    Vector *homols;
    sprintf(qStr,
               "select   grm.member_stable_id,"
               "         gd.name,"
               "         grm.chromosome,"
               "         grm.chrom_start,"
               "         grm.chrom_end"
               " from    gene_relationship_member grm,"
               "         genome_db gd"
               " where   grm.gene_relationship_id = " IDFMTSTR 
               " and     grm.genome_db_id = gd.genome_db_id "
               " and NOT (grm.member_stable_id = '%s')", relationshipIds[i], gene);

    homols = HomologyAdaptor_getHomologues(ha, qStr);
    Vector_append(genes,homols);
    Vector_free(homols);
  }

  free(relationshipIds);
  free(species);

  return genes;
}

Vector **HomologyAdaptor_fetchHomologuesByChrStartInSpecies(HomologyAdaptor *ha, 
                   char *sp, char *chr, int start, char *hSp, int num) {
  char qStr[1024];
  Vector **genesPair = NULL;
  int i;
  char *hSpecies;
  char *species;
  int nRelationship;
  IDType *relationshipIds;

  species = StrUtil_copyString(&species,sp,0);
  species = StrUtil_strReplChr(species,'_',' ');

  hSpecies = StrUtil_copyString(&hSpecies,hSp,0);
  hSpecies = StrUtil_strReplChr(hSpecies,'_',' ');

  sprintf(qStr, 
          "select  grm.gene_relationship_id"
          " from   gene_relationship_member grm,"
          "          genome_db gd"
          "   where  gd.genome_db_id = grm.genome_db_id"
          "   and    gd.name = '%s'"
          "   and    grm.chromosome = '%s'"
          "   and    grm.chrom_start >= %d"
          "   group by grm.gene_relationship_id"
          "   order by grm.chrom_start"
          "   limit %d", species, chr, start, num);

  nRelationship = HomologyAdaptor_getRelationships(ha,qStr,&relationshipIds);

  if ((genesPair = calloc(2,sizeof(Vector *))) == NULL) {
    fprintf(stderr,"Error: Failed allocating genesPair\n");
  } else {
    genesPair[0] = Vector_new();
    genesPair[1] = Vector_new();

    for (i=0;i<nRelationship;i++) {
      Vector *homols;

      homols = HomologyAdaptor_fetchHomologuesBySpeciesRelationshipId(ha, species, relationshipIds[i]);
      Vector_append(genesPair[0],homols);
      Vector_free(homols);

      homols = HomologyAdaptor_fetchHomologuesBySpeciesRelationshipId(ha, hSpecies, relationshipIds[i]);
      Vector_append(genesPair[1],homols);
      Vector_free(homols);
    }
  }

  free(species);
  free(hSpecies);

  return genesPair;
}                               

Vector *HomologyAdaptor_listStableIdsFromSpecies(HomologyAdaptor *ha, char *sp)  {
  StatementHandle *sth;
  ResultRow *row;
  Vector *genes;
  char qStr[1024];
  char *species;

  species = StrUtil_copyString(&species,sp,0);
  species = StrUtil_strReplChr(species,'_',' ');

  sprintf(qStr,
          "select  distinct grm.member_stable_id "
          " from    gene_relationship_member grm,"
          "         genome_db gd "
          " where   gd.genome_db_id = grm.genome_db_id "
          " and     gd.name = '%s'", species);


  sth = ha->prepare((BaseAdaptor *)ha, qStr, strlen(qStr));
  sth->execute(sth);

  genes = Vector_new();

  while ((row = sth->fetchRow(sth))) {
    char *tmpStr;
    Vector_addElement(genes,StrUtil_copyString(&tmpStr, row->getStringAt(row,0),0));
  }
  sth->finish(sth);

  free(species);

  return genes;
}                               


Vector *HomologyAdaptor_fetchHomologuesBySpeciesRelationshipId(HomologyAdaptor *ha, 
                   char *hSpecies, IDType internalId){
  char qStr[1024];

  sprintf(qStr,
           "select  grm.member_stable_id,"
           "        gd.name,"
           "        grm.chromosome,"
           "        grm.chrom_start,"
           "        grm.chrom_end"
           " from    gene_relationship_member grm,"
           "        genome_db gd "
           " where   gd.genome_db_id = grm.genome_db_id "
           " and     gd.name = '%s' "
           " and     grm.gene_relationship_id = " IDFMTSTR,
         hSpecies,internalId);

  return HomologyAdaptor_getHomologues(ha, qStr);
}

Vector *HomologyAdaptor_getHomologues(HomologyAdaptor *ha, char *qStr) {
  StatementHandle *sth;
  ResultRow *row;
  Vector *genes;

  sth = ha->prepare((BaseAdaptor *)ha, qStr, strlen(qStr));
  sth->execute(sth);

  genes = Vector_new();
  while ((row = sth->fetchRow(sth))) {
    Homology *homol = Homology_new();
    Homology_setSpecies(homol, row->getStringAt(row,1));
    Homology_setStableId(homol, row->getStringAt(row,0));
    Homology_setChromosome(homol, row->getStringAt(row,2));
    Homology_setChrStart(homol, row->getIntAt(row,3));
    Homology_setChrEnd(homol, row->getIntAt(row,4));
        
    Vector_addElement(genes,homol);
  }

  sth->finish(sth);

  return genes;
}


int HomologyAdaptor_getRelationships(HomologyAdaptor *ha, char *qStr, IDType **idsP) {
  StatementHandle *sth;
  ResultRow *row;
  int nAlloced = 2;
  int nId = 0;
  int ok = 1;

  sth = ha->prepare((BaseAdaptor *)ha, qStr, strlen(qStr));
  sth->execute(sth);

  if ((*idsP = (IDType *)calloc(nAlloced,sizeof(IDType))) == NULL) {
    fprintf(stderr,"Error: Failed allocating idsP\n");
    ok = 0;
  }

  if (ok) {
    while ((row = sth->fetchRow(sth))) {
      if (nId == nAlloced) {
        nAlloced = (nAlloced *2);
        if ((*idsP = (IDType *)realloc(*idsP,nAlloced*sizeof(IDType))) == NULL) {
          fprintf(stderr,"Error: Failed reallocating idsP\n");
          ok = 0;
          break;
        }
      }

      (*idsP)[nId++] = row->getLongLongAt(row,0); 
    }
  }

  if (ok) {
    sth->finish(sth);
  }

  return nId;
}


IDType HomologyAdaptor_getRelationship(HomologyAdaptor *ha, char *qStr) {
  StatementHandle *sth;
  ResultRow *row;
  IDType relId;
  

  sth = ha->prepare((BaseAdaptor *)ha, qStr, strlen(qStr));
  sth->execute(sth);

  if ((row = sth->fetchRow(sth))) {
    relId = row->getLongLongAt(row,0);
  } else {
    fprintf(stderr, "Error: No relationship found by %s\n",qStr);
    relId = 0;
  }

  sth->finish(sth);

  return relId;
}
