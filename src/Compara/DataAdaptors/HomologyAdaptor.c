#include "HomologyAdaptor.h"


HomologyAdaptor *HomologyAdaptor_new() {
}

package Bio::EnsEMBL::Compara::DBSQL::HomologyAdaptor;
use vars qw(@ISA);
use strict;
use Bio::Root::Object;
use DBI;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Compara::Homology;
@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

=head2 fetch_homologues_of_gene_in_species

 Title   : fetch_homologues_of_gene_in_species
 Usage   : $db->fetch_homologues_of_gene_in_species('Homo_sapiens','ENSG00000218116','Mus_musculus')
 Function: finds homologues, in a certain species, of a given gene
 Example :
 Returns : array of homology objects (Bio::EnsEMBL::Compara::Homology)
 Args    : species gene is from, gene stable name, species to find homology in

=cut

Vector *HomologyAdaptor_fetchHomologuesOfGeneInSpecies(HomologyAdaptor *ha, 
                          char *species, char *gene, char *hSpecies) {
  char qStr[1024];
  Vector *genes;
  IDType *relationshipIds;
  int nRelationship;
  int i;

  $hspecies =~ tr/_/ /;
  $species  =~ tr/_/ /;

  sprintf(qStr,
           "select grm.gene_relationship_id "
           " from   gene_relationship_member grm, "
           "        genome_db gd "
           " where  gd.genome_db_id = grm.genome_db_id "
           " and    gd.name = '%s' "
           " and    grm.member_stable_id = '%s' "
           " group by grm.gene_relationship_id", species, gene);

  nRelationship = HomologyAdaptor_getRelationships(qStr,&relationshipIds);

  genes = Vector_new();
  for (i=0;i<nRelationship;i++) {
    Vector *homols = HomologyAdaptor_fetchHomologuesBySpeciesRelationshipId(ha,hSpecies,relationshipIds[i]);
    Vector_append(genes, homols);
    Vector_free(homols,NULL);
  }

  free(relationshipIds);

  return genes;
}


=head2 fetch_homologues_of_gene

 Title   : fetch_homologues_of_gene
 Usage   : $db->fetch_homologues_of_gene('Homo_sapiens','ENSG00000218116')
 Function: finds homologues of a given gene
 Example :
 Returns : an array of homology objects
 Args    : species gene is from, gene stable name

=cut

Vector *HomologyAdaptor_fetchHomologuesOfGene(HomologyAdaptor *ha, char *species, char *gene) {
  char qStr[1024];

    $species  =~ tr/_/ /;

  sprintf(qStr,
            "select grm.gene_relationship_id "
            " from   gene_relationship_member grm, "
            "        genome_db gd "
            " where  gd.genome_db_id = grm.genome_db_id "
            " and    gd.name = '%s' "
            " and    grm.member_stable_id = '%s' "
            " group by grm.gene_relationship_id", species, gene);

  nRelationship = HomologyAdaptor_getRelationships(qStr,&relationshipIds);

  genes = Vector_new();
  for (i=0;i<nRelationship;i++) {
    Vector *homols;
    sprintf(qStr,
               "select   grm.member_stable_id,
               "         gd.name,
               "         grm.chromosome,
               "         grm.chrom_start,
               "         grm.chrom_end
               " from    gene_relationship_member grm,
               "         genome_db gd
               " where   grm.gene_relationship_id = " IDFMTSTR 
               " and     grm.genome_db_id = gd.genome_db_id 
               " and NOT (grm.member_stable_id = '%s')", relationshipIds[i], gene);

    homols = HomologyAdaptor_getHomologues(ha, qStr);
    Vector_append(genes,homols);
    Vector_free(homols, NULL);
  }

  free(relationshipIds);

  return genes;
}


=head2 fetch_homologues_by_chr_start_in_species

 Title   : fetch_homologues_by_chr_start_in_species
 Usage   : $db->fetch_homologues_by_chr_start_in_species(   'Homo_sapiens',
                                                            'X',
                                                            1000,
                                                            'Mus_musculus', 
                                                            10)
        Will return 10 human genes in order from 1000bp on chrX along with 10 homologues from mouse.  If no homologue, homologue will be "no known homologue"
 Function: finds a list of homologues with a given species
 Example :
 Returns : a hash of species names against arrays of homology objects
 Args    : species from, chr, start, second spp, number of homologues to fetch

=cut

Vector **HomologyAdaptor_fetchHomologuesByChrStartInSpecies(HomologyAdaptor *ha, 
                   char *species, char *chr, int start, char *hspecies, int num) {
  char qStr[1024];
  Vector **genesPair;
  int i;

    $hspecies =~ tr/_/ /;
    $species =~ tr/_/ /;

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

  nRelationship = HomologyAdaptor_getRelationships(qStr,&relationshipIds);

  if ((genesPair = calloc(2,sizeof(Vector *))) == NULL) {
    fprintf(stderr,"Error: Failed allocating genesPair\n");
    exit(1);
  }

  genesPair[0] = Vector_new();
  genesPair[1] = Vector_new();

  for (i=0;i<nRelationship;i++) {
    Vector *homols;

    homols = HomologyAdaptor_fetchHomologuesBySpeciesRelationshipId(ha, species, relationshipIds[i]);
    Vector_append(genesPair[0],homols);
    Vector_free(homols, NULL);

    homols = HomologyAdaptor_fetchHomologuesBySpeciesRelationshipId(ha, hSpecies, relationshipIds[i]);
    Vector_append(genesPair[1],homols);
    Vector_free(homols, NULL);
  }

  return genesPair;
}                               



=head2 list_stable_ids_from_species

 Title   : list_stable_ids_from_species
 Usage   : $db->list_stable_ids_from_species('Homo_sapiens')
 Function: Find all the stable ids in the gene_relationship_member table 
           from a specified species 
 Example :
 Returns : array from transcript stable ids
 Args    : species

=cut

sub HomologyAdaptor_listStableIdsFromSpecies(HomologyAdaptor *ha, char *species)  {
  StatementHandle *sth;
  ResultRow *row;
  Vector *genes;


    $species =~ tr/_/ /;

  char qStr[1024];

  sprintf(qStr,
          "select  grm.member_stable_id "
          " from    gene_relationship_member grm,"
          "         genome_db gd "
          " where   gd.genome_db_id = grm.genome_db_id "
          " and     gd.name = '%s'", species);


  sth = ((BaseAdaptor *)ha)->prepare(qStr);
  sth->execute(sth);

  genes = Vector_new();

  while ((row = sth->fetchRow(sth))) {
    char *tmpStr;
    Vector_addElement(genes,StrUtil_copyString(&tmpStr, row->getStringAt(row,0),0));
  }
  sth->finish(sth);

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
         hSpecies,internalID);

  return HomologyAdaptor_getHomologues(ha, qStr);
}

Vector *HomologyAdaptor_getHomologues(HomologyAdaptor *ha, char *qStr) {
  StatementHandle *sth;
  ResultRow *row;
  Vector *genes;

  sth = ((BaseAdaptor *)ha)->prepare(qStr);
  sth->execute(sth);

  genes = Vector_new();
  while ((row = sth->fetchRow(sth))) {
    Homology *homol = Homology_new();
    Homology_setSpecies(homol, row->getStringAt(->species($ref->{'name'});
    Homology_setStableId(homol, row->getStringAt(->stable_id($ref->{'member_stable_id'});
    Homology_setChromosome(homol, row->getStringAt(->chromosome($ref->{'chromosome'});
    Homology_setChrStart(homol, row->getIntAt(->chrom_start($ref->{'chrom_start'});
    Homology_setChrEnd(homol, row->getIntAt(->chrom_end($ref->{'chrom_end'});
        
    Vector_addElement(genes,homol);
  }

  sth->finish(sth);

  return genes;
}


int HomologyAdaptor_getRelationships(HomologyAdaptor *ha, char *qStr, IDType **idsP) {
  StatementHandle *sth;
  ResultRow *row;
  Vector *genes;
  int nAlloced = 2;
  int nId = 0;

  sth = ((BaseAdaptor *)ha)->prepare(qStr);
  sth->execute(sth);

  if ((*idsP = (IDType *)calloc(nAlloced,sizeof(IDType))) == NULL) {
    fprintf(stderr,"Error: Failed allocating idsP\n");
    exit(1);
  }


  while ((row = sth->fetchRow(sth))) {
    if (nId == nAlloced) {
      nAlloced = (nAlloced *2);
      if ((*idsP = (IDType *)realloc(*idsP,nAlloced*sizeof(IDType))) == NULL) {
        fprintf(stderr,"Error: Failed reallocating idsP\n");
        exit(1);
      }
    }

    (*idsP)[nId++] = row->getLongLongAt(row,0)); 
  }
  sth->finish(sth);

  return nId;
}


IDType HomologyAdaptor_getRelationship(HomologyAdaptor *ha, char *qStr) {
  StatementHandle *sth;
  ResultRow *row;
  IDType relId;
  

  sth = ((BaseAdaptor *)ha)->prepare(qStr);
  sth->execute(sth);

  if ((row = sth->fetchRow(sth))) {
    relId = row->getLongLongAt(row,0);
  } else {
    fprintf(stderr, "Error: No relationship found by %s\n,qStr);
    relId = 0;
  }

  sth->finish(sth);

  return relId;
}
