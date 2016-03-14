/*
 * Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

/*
=head1 DESCRIPTION

This is a database aware adaptor for the retrieval and storage of gene
objects.

=head1 METHODS
*/

#include "GeneAdaptor.h"
#include "AnalysisAdaptor.h"
#include "AssemblyMapperAdaptor.h"
#include "DBAdaptor.h"
#include "DBEntryAdaptor.h"
#include "StrUtil.h"
#include "BaseAdaptor.h"
#include "MysqlUtil.h"
#include "Exon.h"
#include "IDHash.h"
#include "Transcript.h"
#include "Slice.h"
#include "AssemblyMapper.h"
#include "ChainedAssemblyMapper.h"
#include "CoordSystemAdaptor.h"
#include "SupportingFeatureAdaptor.h"
#include "MetaCoordContainer.h"

#include "ExonAdaptor.h"
#include "TranscriptAdaptor.h"
#include "SliceAdaptor.h"
#include "AttributeAdaptor.h"

#include "StatementHandle.h"
#include "ResultRow.h"

#include "Error.h"

#include "ProcUtil.h"

NameTableType GeneAdaptor_tableNames = {{"gene","g"},
                                        {"xref","x"},
                                        {"external_db","exdb"},
                                        {NULL,NULL}};

NameTableType GeneAdaptor_leftJoins = {{"xref","x.xref_id = g.display_xref_id"},
                                       {"external_db","exdb.external_db_id = x.external_db_id"},
                                       {NULL,NULL}};

GeneAdaptor *GeneAdaptor_new(DBAdaptor *dba) {
  GeneAdaptor *ga;

  if ((ga = (GeneAdaptor *)calloc(1,sizeof(GeneAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for GeneAdaptor\n");
    return NULL;
  }
  BaseFeatureAdaptor_init((BaseFeatureAdaptor *)ga, dba, GENE_ADAPTOR);
  //?? ga->sliceGeneCache = StringHash_new(STRINGHASH_MEDIUM);

  ga->getTables                  = GeneAdaptor_getTables;
  ga->getColumns                 = GeneAdaptor_getColumns;
  ga->store                      = GeneAdaptor_store;
  ga->objectsFromStatementHandle = (BaseAdaptor_ObjectsFromStatementHandleFunc)GeneAdaptor_objectsFromStatementHandle;
  ga->leftJoin                   = GeneAdaptor_leftJoin;

  return ga;
}

/*
# _tables
#  Arg [1]    : none
#  Description: PROTECTED implementation of superclass abstract method.
#               Returns the names, aliases of the tables to use for queries.
#  Returntype : list of listrefs of strings
#  Exceptions : none
#  Caller     : internal
#  Status     : Stable
*/
NameTableType *GeneAdaptor_getTables() {
  return &GeneAdaptor_tableNames;
}

/*
# _columns
#  Arg [1]    : none
#  Example    : none
#  Description: PROTECTED implementation of superclass abstract method.
#               Returns a list of columns to use for queries.
#  Returntype : list of strings
#  Exceptions : none
#  Caller     : internal
#  Status     : Stable
*/

// Perl did some extra stuff around the dates to allow cross rdbm portability, but I'm not bothering with that because it will make setting
// up the array much more painfull

char *Gene_cols[] = {
                "g.gene_id",
                "g.seq_region_id",
                "g.seq_region_start",
                "g.seq_region_end",
                "g.seq_region_strand",
                "g.analysis_id",
                "g.biotype",
                "g.display_xref_id",
                "g.description",
                "g.status",
                "g.source",
                "g.is_current",
                "g.canonical_transcript_id",
                //"g.canonical_annotation",
                "g.stable_id",
                "g.version",
                "UNIX_TIMESTAMP(g.created_date)",
                "UNIX_TIMESTAMP(g.modified_date)",
                "x.display_label",
                "x.dbprimary_acc",
                "x.description",
                "x.version",
                "exdb.db_name",
                "exdb.status",
                "exdb.db_release",
                "exdb.db_display_name",
                "x.info_type",
                "x.info_text",
                NULL };

char **GeneAdaptor_getColumns() {
  return Gene_cols;
}

NameTableType *GeneAdaptor_leftJoin() {
  return &GeneAdaptor_leftJoins;
}

/*
=head2 list_dbIDs

  Example    : @gene_ids = @{$gene_adaptor->list_dbIDs()};
  Description: Gets an array of internal ids for all genes in the current db
  Arg[1]     : <optional> int. not 0 for the ids to be sorted by the seq_region.
  Returntype : Listref of Ints
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/
Vector *GeneAdaptor_listDbIDs(GeneAdaptor *ga, int ordered) {

// NIY: Shouldn't really do this direct BaseAdaptor call but ...
  return BaseAdaptor_listDbIDs((BaseAdaptor *)ga, "gene", NULL, ordered);
}

/*
=head2 list_stable_ids

  Example    : @stable_gene_ids = @{$gene_adaptor->list_stable_ids()};
  Description: Gets an listref of stable ids for all genes in the current db
  Returntype : reference to a list of strings
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/
Vector *GeneAdaptor_listStableIDs(GeneAdaptor *ga) {

// NIY: Shouldn't really do this direct BaseAdaptor call but ...
  return BaseAdaptor_listDbIDs((BaseAdaptor *)ga, "gene", "stable_id", 0);
}

Vector *GeneAdaptor_listSeqRegionIds(GeneAdaptor *ga) {

// NIY: Shouldn't really do this direct BaseFeatureAdaptor call but ...
  return BaseFeatureAdaptor_listSeqRegionIds((BaseFeatureAdaptor *)ga, "gene");
}

/*
=head2 fetch_by_display_label

  Arg [1]    : String $label - display label of gene to fetch
  Example    : my $gene = $geneAdaptor->fetch_by_display_label("BRCA2");
  Description: Returns the gene which has the given display label or undef if
               there is none. If there are more than 1, the gene on the 
               reference slice is reported or if none are on the reference,
               the first one is reported.
  Returntype : Bio::EnsEMBL::Gene
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/
Gene *GeneAdaptor_fetchByDisplayLabel(GeneAdaptor *ga, char *label) {
  char constraint[1024];

  sprintf(constraint, "x.display_label = '%s' AND g.is_current = 1", label);

  Vector *genes = GeneAdaptor_genericFetch(ga, constraint, NULL, NULL);

  Gene *gene = NULL;

  if (Vector_getNumElement(genes) > 1) {
    int i;
    for (i=0; i<Vector_getNumElement(genes) && gene == NULL; i++) {
      Gene *g = Vector_getElementAt(genes, i);
     
      if (Slice_isReference(Gene_getSlice(g))) {
        gene = g;
      }
    }
    if (gene == NULL) {
      gene = Vector_getElementAt(genes, 0);
    }
  } else if (Vector_getNumElement(genes) == 1) {
    gene = Vector_getElementAt(genes, 0);
  }

  Vector_free(genes);

  return gene;
}

/*
=head2 fetch_all_by_display_label

  Arg [1]    : String $label - display label of genes to fetch
  Example    : my @genes = @{$geneAdaptor->fetch_all_by_display_label("PPP1R2P1")};
  Description: Returns all genes which have the given display label or undef if
               there are none. 
  Returntype : listref of Bio::EnsEMBL::Gene objects
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/
Vector *GeneAdaptor_fetchAllByDisplayLabel(GeneAdaptor *ga, char *label) {
  char constraint[1024];

  sprintf(constraint, "x.display_label = '%s' AND g.is_current = 1", label);

  return GeneAdaptor_genericFetch(ga, constraint, NULL, NULL);
}

/*
=head2 fetch_by_stable_id

  Arg [1]    : String $id 
               The stable ID of the gene to retrieve
  Example    : $gene = $gene_adaptor->fetch_by_stable_id('ENSG00000148944');
  Description: Retrieves a gene object from the database via its stable id.
               The gene will be retrieved in its native coordinate system (i.e.
               in the coordinate system it is stored in the database). It may
               be converted to a different coordinate system through a call to
               transform() or transfer(). If the gene or exon is not found
               undef is returned instead.
  Returntype : Bio::EnsEMBL::Gene or undef
  Exceptions : if we cant get the gene in given coord system
  Caller     : general
  Status     : Stable

=cut
*/
Gene *GeneAdaptor_fetchByStableId(GeneAdaptor *ga, char *stableId) {
  char constraint[1024];

  sprintf(constraint, "g.stable_id = '%s' AND g.is_current = 1", stableId);

  Vector *genes = GeneAdaptor_genericFetch(ga, constraint, NULL, NULL);
  Gene *gene = Vector_getElementAt(genes, 0);
  Vector_free(genes);

  return gene;
}

/*
=head2 fetch_all_by_biotype 

  Arg [1]    : String $biotype 
               listref of $biotypes
               The biotype of the gene to retrieve. You can have as an argument a reference
               to a list of biotypes
  Example    : $gene = $gene_adaptor->fetch_all_by_biotype('protein_coding'); 
               $gene = $gene_adaptor->fetch_all_by_biotypes(['protein_coding', 'sRNA', 'miRNA']);
  Description: Retrieves an array reference of gene objects from the database via its biotype or biotypes.
               The genes will be retrieved in its native coordinate system (i.e.
               in the coordinate system it is stored in the database). It may
               be converted to a different coordinate system through a call to
               transform() or transfer(). If the gene or exon is not found
               undef is returned instead.
  Returntype  : listref of Bio::EnsEMBL::Gene
  Exceptions : if we cant get the gene in given coord system
  Caller     : general
  Status     : Stable

=cut
*/
Vector *GeneAdaptor_fetchAllByBiotype(GeneAdaptor *ga, Vector *biotypes) {
  char constraint[1024];

  GeneAdaptor_biotypeConstraint(ga, biotypes, constraint);
  
  return GeneAdaptor_genericFetch(ga, constraint, NULL, NULL);
}


/*
=head2 biotype_constraint 

  Arg [1]    : String $biotype 
               listref of $biotypes
               The biotype of the gene to retrieve. You can have as an argument a reference
               to a list of biotypes
  Description: Used internally to generate a SQL constraint to restrict a gene query by biotype
  Returntype  : String
  Exceptions : If biotype is not supplied
  Caller     : general
  Status     : Stable

=cut
*/
void GeneAdaptor_biotypeConstraint(GeneAdaptor *ga, Vector *biotypes, char *constraint) {
  
  if (biotypes == NULL || Vector_getNumElement(biotypes) == 0) {
    fprintf(stderr,"list of biotypes expected\n");
    exit(1);
  }

  if (Vector_getNumElement(biotypes) > 1) {
    strcpy(constraint, "g.biotype IN (");

    int i;
    for (i=0;i<Vector_getNumElement(biotypes); i++) {
      char *biotype = Vector_getElementAt(biotypes, i);
      
      if (i>0) {
        strcat(constraint, ", ");
      }
      sprintf(constraint, "%s'%s'", constraint, biotype);
    }
    strcat(constraint, ") and g.is_current = 1");

  } else { // just one
    char *biotype = Vector_getElementAt(biotypes, 0);

    sprintf(constraint, "g.biotype = '%s' and g.is_current = 1", biotype);
  }

  return;
}

/*
=head2 count_all_by_biotype 

  Arg [1]     : String $biotype 
                listref of $biotypes
                The biotype of the gene to retrieve. You can have as an argument a reference
                to a list of biotypes
  Example     : $cnt = $gene_adaptor->count_all_by_biotype('protein_coding'); 
                $cnt = $gene_adaptor->count_all_by_biotypes(['protein_coding', 'sRNA', 'miRNA']);
  Description : Retrieves count of gene objects from the database via its biotype or biotypes.
  Returntype  : integer
  Caller      : general
  Status      : Stable

=cut
*/
int GeneAdaptor_countAllByBiotype(GeneAdaptor *ga, Vector *biotypes) {
  char constraint[1024];

  GeneAdaptor_biotypeConstraint(ga, biotypes, constraint);
  
  return GeneAdaptor_genericCount(ga, constraint);
}

Vector *GeneAdaptor_fetchAll(GeneAdaptor *ga) {
  char constraint[1024];

// NIY: Maybe just a constant string
  sprintf(constraint, "g.biotype != 'LRG_gene' and g.is_current = 1");

  return GeneAdaptor_genericFetch(ga, constraint, NULL, NULL);
}

/*
=head2 fetch_all_versions_by_stable_id 

  Arg [1]     : String $stable_id 
                The stable ID of the gene to retrieve
  Example     : $gene = $gene_adaptor->fetch_all_versions_by_stable_id
                  ('ENSG00000148944');
  Description : Similar to fetch_by_stable_id, but retrieves all versions of a
                gene stored in the database.
  Returntype  : listref of Bio::EnsEMBL::Gene
  Exceptions  : if we cant get the gene in given coord system
  Caller      : general
  Status      : At Risk

=cut
*/
Vector *GeneAdaptor_fetchAllVersionsByStableId(GeneAdaptor *ga, char *stableId) {
  char constraint[1024];

  sprintf(constraint, "g.stable_id = '%s'", stableId);

  return GeneAdaptor_genericFetch(ga, constraint, NULL, NULL);
}

/*
=head2 fetch_by_exon_stable_id

  Arg [1]    : String $id
               The stable id of an exon of the gene to retrieve
  Example    : $gene = $gene_adptr->fetch_by_exon_stable_id('ENSE00000148944');
  Description: Retrieves a gene object from the database via an exon stable id.
               The gene will be retrieved in its native coordinate system (i.e.
               in the coordinate system it is stored in the database). It may
               be converted to a different coordinate system through a call to
               transform() or transfer(). If the gene or exon is not found
               undef is returned instead.
  Returntype : Bio::EnsEMBL::Gene or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/
// Note perl had a version arg, but it didn't seem to be used????
Gene *GeneAdaptor_fetchByExonStableId(GeneAdaptor *ga, char *stableId) {
  char qStr[1024];

  sprintf(qStr,
      "SELECT t.gene_id "
        "FROM transcript as t, "
             "exon_transcript as et, "
             "exon as e "
       "WHERE t.transcript_id = et.transcript_id "
         "AND et.exon_id = e.exon_id "
         "AND e.stable_id = '%s' "
         "AND e.is_current = 1", stableId);

  StatementHandle *sth = ga->prepare((BaseAdaptor *)ga,qStr,strlen(qStr));
  sth->execute(sth);

  if (sth->numRows(sth) != 1) {
    fprintf(stderr, "Failed fetching gene using exon stable id %s - returning NULL\n",  stableId);
    return NULL;
  }

  ResultRow *row = sth->fetchRow(sth);

  IDType dbID = row->getLongLongAt(row, 0);
  sth->finish(sth);

  Gene *gene = (Gene *)GeneAdaptor_fetchByDbID(ga, dbID);

  return gene;
}

/* Don't bother
=head2 fetch_all_by_domain

  Arg [1]    : String $domain
               The domain to fetch genes from
  Example    : my @genes = @{ $gene_adaptor->fetch_all_by_domain($domain) };
  Description: Retrieves a listref of genes whose translation contain interpro
               domain $domain. The genes are returned in their native coord
               system (i.e. the coord_system they are stored in). If the coord
               system needs to be changed, then tranform or transfer should be
               called on the individual objects returned.
  Returntype : list of Bio::EnsEMBL::Genes
  Exceptions : none
  Caller     : domainview
  Status     : Stable

=cut

sub fetch_all_by_domain {
  my ($self, $domain) = @_;

  throw("domain argument is required") unless ($domain);

  my $sth = $self->prepare(
  qq(
  SELECT    tr.gene_id
  FROM      interpro i,
            protein_feature pf,
            transcript tr,
            translation tl,
            seq_region sr,
            coord_system cs
  WHERE     cs.species_id = ?
    AND     cs.coord_system_id = sr.coord_system_id
    AND     sr.seq_region_id = tr.seq_region_id
    AND     tr.is_current = 1
    AND     tr.transcript_id = tl.transcript_id
    AND     tl.translation_id = pf.translation_id
    AND     pf.hit_name = i.id
    AND     i.interpro_ac = ?
  GROUP BY  tr.gene_id));

  $sth->bind_param(1, $self->species_id(), SQL_VARCHAR);
  $sth->bind_param(2, $domain,             SQL_VARCHAR);

  $sth->execute();

  my @array = @{$sth->fetchall_arrayref()};
  $sth->finish();

  my @gene_ids = map { $_->[0] } @array;

  return $self->fetch_all_by_dbID_list(\@gene_ids);
} ## end sub fetch_all_by_domain
*/

/*
=head2 fetch_all_by_Slice_and_external_dbname_link

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               The slice to fetch genes on.
  Arg [2]    : (optional) string $logic_name
               the logic name of the type of features to obtain
  Arg [3]    : (optional) boolean $load_transcripts
               if true, transcripts will be loaded immediately
               rather than lazy loaded later.
  Arg [4]    : Name of the external database
  Example    : @genes = @{
                 $ga->fetch_all_by_Slice_and_external_dbname_link(
                                          $slice, undef, undef, "HUGO" ) };
  Description: Overrides superclass method to optionally load
               transcripts immediately rather than lazy-loading them
               later.  This is more efficient when there are a lot
               of genes whose transcripts are going to be used. The
               genes are then filtered to return only those with
               external database links of the type specified
  Returntype : reference to list of genes
  Exceptions : thrown if exon cannot be placed on transcript slice
  Caller     : 
  Status     : Stable

=cut
*/
/* NIY
Vector *GeneAdaptor_fetchAllBySliceAndExternalDbnameLink(GeneAdaptor *ga, Slice *slice, char *logicName, int loadTranscripts, char *dbName) {

  // Get the external_db_id(s) from the name.
  my $sth = $self->prepare("SELECT external_db_id FROM external_db WHERE db_name = ?");

  $sth->bind_param(1, $db_name, SQL_VARCHAR);
  $sth->execute();

  my $external_db_id;
  $sth->bind_columns(\$external_db_id);

  my @external_db_ids;
  while ($sth->fetch()) {
    push(@external_db_ids, $external_db_id);
  }

  if (scalar(@external_db_ids) == 0) {
    warn sprintf("Could not find external database " . "'%s' in the external_db table\n" . "Available are:\n", $db_name);

    $sth = $self->prepare("SELECT DISTINCT db_name FROM external_db");

    $sth->execute();
    $sth->bind_columns(\$external_db_id);

    while ($sth->fetch()) {
      warn "\t$external_db_id\n";
    }
    return [];
  }

  // Get the gene_ids for those with links.
  my $dbe_adaptor = $self->db()->get_DBEntryAdaptor();

  my %linked_genes;
  foreach $external_db_id (@external_db_ids) {
    my @linked_genes = $dbe_adaptor->list_gene_ids_by_external_db_id($external_db_id);

    foreach my $gene_id (@linked_genes) {
      $linked_genes{$gene_id} = 1;
    }
  }

  // Get all the genes on the slice.
  my $genes = $self->SUPER::fetch_all_by_Slice_constraint($slice, 'g.is_current = 1', $logic_name);

  // Create a list of those that are in the gene_ids list.
  my @genes_passed;
  foreach my $gene (@$genes) {
    if (exists($linked_genes{$gene->dbID()})) {
      push(@genes_passed, $gene);
    }
  }

  // Return the list of those that passed.
  return \@genes_passed;
}
*/

/*
=head2 fetch_all_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               The slice to fetch genes on.
  Arg [2]    : (optional) string $logic_name
               the logic name of the type of features to obtain
  Arg [3]    : (optional) boolean $load_transcripts
               if true, transcripts will be loaded immediately rather than
               lazy loaded later.
  Arg [4]    : (optional) string $source
               the source name of the features to obtain.
  Arg [5]    : (optional) string biotype
                the biotype of the features to obtain.
  Example    : @genes = @{$gene_adaptor->fetch_all_by_Slice()};
  Description: Overrides superclass method to optionally load transcripts
               immediately rather than lazy-loading them later.  This
               is more efficient when there are a lot of genes whose
               transcripts are going to be used.
  Returntype : reference to list of genes 
  Exceptions : thrown if exon cannot be placed on transcript slice
  Caller     : Slice::get_all_Genes
  Status     : Stable

=cut
*/

Vector *GeneAdaptor_fetchAllBySlice(GeneAdaptor *ga, Slice *slice, char *logicName, int loadTranscripts, char *source, char *biotype) {
  char constraint[1024];

  strcpy(constraint, "g.is_current = 1");

  if (source != NULL) {
    sprintf(constraint,"%s and g.source = '%s'", constraint, source);
  }
  if (biotype != NULL) {
    sprintf(constraint,"%s and g.biotype = '%s'", constraint, biotype);
  }

  // Perl did a direct SUPER::fetch_all_by_Slice_constraint - not sure why???

  Vector *genes = GeneAdaptor_fetchAllBySliceConstraint(ga, slice, constraint, logicName);

  // If there are less than two genes, still do lazy-loading.
// SMJS Tweaked so never does lazy loading
  if (!loadTranscripts || Vector_getNumElement(genes) < 1) {
    return genes;
  }


  // Preload all of the transcripts now, instead of lazy loading later,
  // faster than one query per transcript.

  // First check if transcripts are already preloaded.
  // FIXME: Should check all transcripts.
  Gene *firstGene = Vector_getElementAt(genes, 0);
  // Deliberate direct reach into gene struct to avoid any lazy loading
  // Not sure if I should set transcripts to null when there aren't any loaded - maybe, but it makes everything harder
  if (Vector_getNumElement(firstGene->transcripts) != 0) { // any transcripts???
// NIY: Is there anything to free??
    return genes;
  }

  // Get extent of region spanned by transcripts.
  long minStart =  2000000000;
  long maxEnd   = -2000000000;

  int i;
  for (i=0; i<Vector_getNumElement(genes); i++) {
    Gene *g  = Vector_getElementAt(genes, i);
    if (Gene_getSeqRegionStart((SeqFeature*)g) < minStart) {
      minStart = Gene_getSeqRegionStart((SeqFeature*)g);
    }
    if (Gene_getSeqRegionEnd((SeqFeature*)g) > maxEnd) {
      maxEnd = Gene_getSeqRegionEnd((SeqFeature*)g);
    }
  }

  Slice *extSlice;

  if (minStart >= Slice_getStart(slice) && maxEnd <= Slice_getEnd(slice)) {
    extSlice = slice;
  } else {
    SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(ga->dba);
    extSlice = SliceAdaptor_fetchByRegion(sa, Slice_getCoordSystemName(slice), Slice_getSeqRegionName(slice), 
                                          minStart, maxEnd, Slice_getStrand(slice), CoordSystem_getVersion(Slice_getCoordSystem(slice)), 0);
  }

  // Associate transcript identifiers with genes.
  IDHash *gHash = IDHash_new(IDHASH_MEDIUM);
  for (i=0; i<Vector_getNumElement(genes); i++) {
    Gene *g  = Vector_getElementAt(genes, i);
    if ( ! IDHash_contains(gHash, Gene_getDbID(g))) {
      IDHash_add(gHash, Gene_getDbID(g), g);
    }
  }

  IDType *uniqueIds = IDHash_getKeys(gHash);

  char *qStr = NULL;
  if ((qStr = (char *)calloc(655500,sizeof(char))) == NULL) {
    fprintf(stderr,"Failed allocating qStr\n");
    return genes;
  }

  strcpy(qStr, "SELECT gene_id, transcript_id FROM   transcript WHERE  gene_id IN (");
  for (i=0; i<IDHash_getNumValues(gHash); i++) {
    if (i!=0) {
      strcat(qStr, ", ");
    }
    sprintf(qStr, "%s"IDFMTSTR, qStr, uniqueIds[i]);
  }
  strcat(qStr,")");

  free(uniqueIds);

  StatementHandle *sth = ga->prepare((BaseAdaptor *)ga,qStr,strlen(qStr));
  sth->execute(sth);

  IDHash *trGHash = IDHash_new(IDHASH_MEDIUM);
  ResultRow *row;
// Extra parentheses to please mac compiler
  while ((row = sth->fetchRow(sth))) {
    IDType gId  = row->getLongLongAt(row,0);
    IDType trId = row->getLongLongAt(row,1);

    IDHash_add(trGHash, trId, IDHash_getValue(gHash, gId));
  }

  IDHash_free(gHash, NULL);

  sth->finish(sth);

  TranscriptAdaptor *ta = DBAdaptor_getTranscriptAdaptor(ga->dba);

  // Note this is a constraint rather than a complete query, but I'm using qStr to save stack space
  uniqueIds = IDHash_getKeys(trGHash);

  qsort(uniqueIds, IDHash_getNumValues(trGHash), sizeof(IDType), idTypeCompFunc);
//sprintf("t.transcript_id IN (%s)", join(',', sort { $a <=> $b } keys(%tr_g_hash))));


  char tmpStr[1024];
  int lenNum;
  int endPoint = sprintf(qStr, "t.transcript_id IN (");
  for (i=0; i<IDHash_getNumValues(trGHash); i++) {
    if (i!=0) {
      qStr[endPoint++] = ',';
      qStr[endPoint++] = ' ';
      //strcat(qStr, ", ");
    }
    //fprintf(stderr, "id %d %s\n", i, tmpStr);
    //sprintf(qStr, "%s"IDFMTSTR, qStr, uniqueIds[i]);
    lenNum = sprintf(tmpStr,IDFMTSTR,uniqueIds[i]);
    memcpy(&(qStr[endPoint]), tmpStr, lenNum);
    endPoint+=lenNum;
  }
  //strcat(qStr,")");
  qStr[endPoint++] = ')';
  qStr[endPoint] = '\0';

  free(uniqueIds);

  Vector *transcripts = TranscriptAdaptor_fetchAllBySlice(ta, extSlice, 1, NULL, qStr /*which is transcript constraint*/);

  // Move transcripts onto gene slice, and add them to genes.
  Vector *exons = Vector_new();
  for (i=0; i<Vector_getNumElement(transcripts); i++) {
    Transcript *tr = Vector_getElementAt(transcripts, i);
    
    if (!IDHash_contains(trGHash, Transcript_getDbID(tr))) continue;

    Transcript *newTr;
    if (slice != extSlice) {
      newTr = Transcript_transfer(tr, slice);
      if (newTr == NULL) {
        fprintf(stderr, "Unexpected. Transcript could not be transferred onto Gene slice.\n");
        exit(1);
      }
    } else {
      newTr = tr;
    }

    // Note perl used tr->dbID rather than new_tr, but I may free that before here
    Gene *g = IDHash_getValue(trGHash, Transcript_getDbID(newTr));
    Gene_addTranscript(g, newTr);

    // Used below for fetching 
    Vector_append(exons, Transcript_getAllExons(newTr));
  }

  IDHash_free(trGHash, NULL);

  SupportingFeatureAdaptor *sfa = DBAdaptor_getSupportingFeatureAdaptor(ga->dba);
  Vector *tmpVec = SupportingFeatureAdaptor_fetchAllByExonList(sfa, exons, slice);
  Vector_free(tmpVec);
  Vector_free(exons);
  free(qStr);

  return genes;
}

/*
=head2 count_all_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               The slice to count genes on.
  Arg [2]    : (optional) biotype(s) string or arrayref of strings 
                the biotype of the features to count.
  Arg [1]    : (optional) string $source
               the source name of the features to count.
  Example    : $cnt = $gene_adaptor->count_all_by_Slice();
  Description: Method to count genes on a given slice, filtering by biotype and source
  Returntype : integer
  Exceptions : thrown if exon cannot be placed on transcript slice
  Status     : Stable
  Caller     : general
=cut
*/
int GeneAdaptor_countAllBySlice(GeneAdaptor *ga, Slice *slice, Vector *biotypes, char *source) {
  char constraint[1024];

  strcpy(constraint, "g.is_current = 1");
  if (source != NULL) {
    sprintf(constraint,"%s and g.source = '%s'", constraint, source);
  }
  if (biotypes != NULL) {
    char biotypeConstraint[1024];
    GeneAdaptor_biotypeConstraint(ga, biotypes, biotypeConstraint);
    sprintf(constraint, "%s AND %s", constraint, biotypeConstraint);
  }

  return GeneAdaptor_countBySliceConstraint(ga, slice, constraint, NULL);
}

/*
=head2 fetch_by_transcript_id

  Arg [1]    : Int $trans_id
               Unique database identifier for the transcript whose gene should
               be retrieved. The gene is returned in its native coord
               system (i.e. the coord_system it is stored in). If the coord
               system needs to be changed, then tranform or transfer should
               be called on the returned object. undef is returned if the
               gene or transcript is not found in the database.
  Example    : $gene = $gene_adaptor->fetch_by_transcript_id(1241);
  Description: Retrieves a gene from the database via the database identifier
               of one of its transcripts.
  Returntype : Bio::EnsEMBL::Gene
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/
Gene *GeneAdaptor_fetchByTranscriptId(GeneAdaptor *ga, IDType transId) {
  char qStr[1024];

  sprintf(qStr,"SELECT tr.gene_id "
               "FROM transcript tr "
               "WHERE tr.transcript_id = "IDFMTSTR, transId);

  StatementHandle *sth = ga->prepare((BaseAdaptor *)ga,qStr,strlen(qStr));
  sth->execute(sth);

  if (sth->numRows(sth) == 0) {
    return NULL;
  }

  ResultRow *row = sth->fetchRow(sth);
  IDType geneId = row->getLongLongAt(row, 0);

  sth->finish(sth);

  Gene *gene = (Gene *)GeneAdaptor_fetchByDbID(ga, geneId);
  return gene;
}

/*
=head2 fetch_by_transcript_stable_id

  Arg [1]    : string $trans_stable_id
               transcript stable ID whose gene should be retrieved
  Example    : my $gene = $gene_adaptor->fetch_by_transcript_stable_id
                 ('ENST0000234');
  Description: Retrieves a gene from the database via the stable ID of one of
               its transcripts
  Returntype : Bio::EnsEMBL::Gene
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/
Gene *GeneAdaptor_fetchByTranscriptStableId(GeneAdaptor *ga, char *transStableId) {
  char qStr[1024];

  sprintf(qStr,"SELECT tr.gene_id "
               "FROM transcript tr "
               "WHERE tr.is_current = 1 AND tr.stable_id = '%s'", transStableId);

  StatementHandle *sth = ga->prepare((BaseAdaptor *)ga,qStr,strlen(qStr));
  sth->execute(sth);

  if (sth->numRows(sth) == 0) {
    return NULL;
  }

  ResultRow *row = sth->fetchRow(sth);
  IDType geneId = row->getLongLongAt(row, 0);

  sth->finish(sth);

  Gene *gene = (Gene *)GeneAdaptor_fetchByDbID(ga, geneId);
  return gene;
}

/*
=head2 fetch_by_translation_stable_id

  Arg [1]    : String $translation_stable_id
               The stable id of a translation of the gene to be obtained
  Example    : my $gene = $gene_adaptor->fetch_by_translation_stable_id
                 ('ENSP00000278194');
  Description: Retrieves a gene via the stable id of one of its translations.
  Returntype : Bio::EnsEMBL::Gene
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/
Gene *GeneAdaptor_fetchByTranslationStableId(GeneAdaptor *ga, char *translationStableId) {
  char qStr[1024];

  sprintf(qStr, "SELECT  tr.gene_id "
                "FROM    transcript tr, translation tl "
                "WHERE   tl.stable_id = '%s' "
                "AND     tr.transcript_id = tl.transcript_id "
                "AND     tr.is_current = 1", translationStableId);

  StatementHandle *sth = ga->prepare((BaseAdaptor *)ga,qStr,strlen(qStr));
  sth->execute(sth);

  if (sth->numRows(sth) == 0) {
    return NULL;
  }

  ResultRow *row = sth->fetchRow(sth);
  IDType geneId = row->getLongLongAt(row, 0);

  sth->finish(sth);

  Gene *gene = (Gene *)GeneAdaptor_fetchByDbID(ga, geneId);
  return gene;
}

/*
=head2 fetch_all_by_external_name

  Arg [1]    : String $external_name
               The external identifier for the gene to be obtained
  Arg [2]    : (optional) String $external_db_name
               The name of the external database from which the
               identifier originates.
  Arg [3]    : Boolean override. Force SQL regex matching for users
               who really do want to find all 'NM%'
  Example    : @genes = @{$gene_adaptor->fetch_all_by_external_name('BRCA2')}
               @many_genes = @{$gene_adaptor->fetch_all_by_external_name('BRCA%')}
  Description: Retrieves a list of genes with an external database
               identifier $external_name. The genes returned are in
               their native coordinate system, i.e. in the coordinate
               system they are stored in the database in.  If another
               coordinate system is required then the Gene::transfer or
               Gene::transform method can be used.
               SQL wildcards % and _ are supported in the $external_name,
               but their use is somewhat restricted for performance reasons.
               Users that really do want % and _ in the first three characters
               should use argument 3 to prevent optimisations
  Returntype : listref of Bio::EnsEMBL::Gene
  Exceptions : none
  Caller     : goview, general
  Status     : Stable

=cut
*/
/* NIY
sub fetch_all_by_external_name {
  my ($self, $external_name, $external_db_name, $override) = @_;

  my $entryAdaptor = $self->db->get_DBEntryAdaptor();

  my @ids = $entryAdaptor->list_gene_ids_by_extids($external_name, $external_db_name, $override);

  my %genes_by_dbIDs =
    map { $_->dbID(), $_ } @{$self->fetch_all_by_dbID_list(\@ids)};

  my @result = map { $genes_by_dbIDs{$_} } @ids;

  return \@result;
}
*/

/* Don't bother
=head2 fetch_all_by_GOTerm

  Arg [1]   : Bio::EnsEMBL::OntologyTerm
              The GO term for which genes should be fetched.

  Example:  @genes = @{
              $gene_adaptor->fetch_all_by_GOTerm(
                $go_adaptor->fetch_by_accession('GO:0030326') ) };

  Description   : Retrieves a list of genes that are associated with
                  the given GO term, or with any of its descendent
                  GO terms.  The genes returned are in their native
                  coordinate system, i.e. in the coordinate system
                  in which they are stored in the database.  If
                  another coordinate system is required then the
                  Gene::transfer or Gene::transform method can be
                  used.

  Return type   : listref of Bio::EnsEMBL::Gene
  Exceptions    : Throws of argument is not a GO term
  Caller        : general
  Status        : Stable

=cut
sub fetch_all_by_GOTerm {
  my ($self, $term) = @_;

  assert_ref($term, 'Bio::EnsEMBL::OntologyTerm');
  if ($term->ontology() ne 'GO') {
    throw('Argument is not a GO term');
  }

  my $entryAdaptor = $self->db->get_DBEntryAdaptor();

  my %unique_dbIDs;
  foreach my $accession (map { $_->accession() } ($term, @{$term->descendants()})) {
    my @ids = $entryAdaptor->list_gene_ids_by_extids($accession, 'GO');
    foreach my $dbID (@ids) { $unique_dbIDs{$dbID} = 1 }
  }

  my @result = @{$self->fetch_all_by_dbID_list([sort { $a <=> $b } keys(%unique_dbIDs)])};

  return \@result;
}

=head2 fetch_all_by_GOTerm_accession

  Arg [1]   : String
              The GO term accession for which genes should be
              fetched.

  Example   :

    @genes =
      @{ $gene_adaptor->fetch_all_by_GOTerm_accession(
        'GO:0030326') };

  Description   : Retrieves a list of genes that are associated with
                  the given GO term, or with any of its descendent
                  GO terms.  The genes returned are in their native
                  coordinate system, i.e. in the coordinate system
                  in which they are stored in the database.  If
                  another coordinate system is required then the
                  Gene::transfer or Gene::transform method can be
                  used.

  Return type   : listref of Bio::EnsEMBL::Gene
  Exceptions    : Throws of argument is not a GO term accession
  Caller        : general
  Status        : Stable

=cut

sub fetch_all_by_GOTerm_accession {
  my ($self, $accession) = @_;

  if ($accession !~ /^GO:/) {
    throw('Argument is not a GO term accession');
  }

  my $goAdaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi', 'Ontology', 'OntologyTerm');

  my $term = $goAdaptor->fetch_by_accession($accession);

  return $self->fetch_all_by_GOTerm($term);
}
*/

/*
=head2 fetch_all_alt_alleles

  Arg [1]    : Bio::EnsEMBL::Gene $gene
               The gene to fetch alternative alleles for
  Example    : my @alt_genes = @{ $gene_adaptor->fetch_all_alt_alleles($gene) };
               foreach my $alt_gene (@alt_genes) {
                 print "Alternate allele: " . $alt_gene->stable_id() . "\n" ;
               }
  Description: Retrieves genes which are alternate alleles to a provided gene.
               Alternate alleles in Ensembl are genes which are similar and are
               on an alternative haplotype of the same region. There are not 
               currently very many of these. This method will return a 
               reference to an empty list if no alternative alleles are found.
  Returntype : listref of Bio::EnsEMBL::Genes
  Exceptions : throw if incorrect arg provided
               warning if gene arg does not have dbID
  Caller     : Gene::get_all_alt_alleles
  Status     : Stable

=cut
*/
Vector *GeneAdaptor_fetchAllAltAlleles(GeneAdaptor *ga, Gene *gene) {
  if (gene == NULL || gene->objectType != CLASS_GENE) {
    fprintf(stderr, "Gene type argument is required\n");
    exit(1);
  }

  IDType geneId = Gene_getDbID(gene);

  if (!geneId) {
    fprintf(stderr, "Warning: Cannot retrieve alternate alleles for gene without dbID\n");
    return Vector_new();
  }
 
  char qStr[1024];
  sprintf(qStr,"SELECT aa1.gene_id "
               "FROM   alt_allele aa1, alt_allele aa2 "
               "WHERE  aa1.alt_allele_id = aa2.alt_allele_id "
               "AND    aa2.gene_id = "IDFMTSTR"AND    aa1.gene_id <> "IDFMTSTR, 
               geneId, geneId);

  StatementHandle *sth = ga->prepare((BaseAdaptor *)ga,qStr,strlen(qStr));
  sth->execute(sth);


  Vector *altIds = Vector_new();
  ResultRow *row;
// Extra parentheses to please mac compiler
  while ((row = sth->fetchRow(sth))) {
    IDType id = row->getLongLongAt(row, 0);
    IDType *idP;

    if ((idP = calloc(1,sizeof(IDType))) == NULL) {
      fprintf(stderr, "Failed allocating space for a id\n");
      exit(1);
    }

    *idP = id;
    Vector_addElement(altIds, idP);
  }
  sth->finish(sth);

  Vector *out;
  if (Vector_getNumElement(altIds) > 0) {
    out = GeneAdaptor_fetchAllByDbIDList(ga, altIds, NULL);
  } else {
    out = Vector_new();
  }

  Vector_setFreeFunc(altIds, free);
  Vector_free(altIds);

  return out;
}

int GeneAdaptor_isRef(GeneAdaptor *ga, IDType geneId) {
  // easier to find if it is not an alt_Allele do this and then negate it.
  char qStr[1024];
  sprintf(qStr,"SELECT count(1) from alt_allele where gene_id = "IDFMTSTR" and is_ref = 0", geneId);
  StatementHandle *sth = ga->prepare((BaseAdaptor *)ga,qStr,strlen(qStr));
  sth->execute(sth);

  if (sth->numRows(sth) == 0) {
    return 1;
  }

  ResultRow *row = sth->fetchRow(sth);
  int isNotRef = row->getIntAt(row, 0);

  if (isNotRef > 0) {
    return 0;
  }

  return 1;
}

/*
=head2 store_alt_alleles


  Arg [1]    : reference to list of Bio::EnsEMBL::Genes $genes
  Example    : $gene_adaptor->store_alt_alleles([$gene1, $gene2, $gene3]);
  Description: This method creates a group of alternative alleles (i.e. locus)
               from a set of genes. The genes should be genes from alternate
               haplotypes which are similar. The genes must already be stored
               in this database. 
  Returntype : int alt_allele_id or undef if no alt_alleles were stored
  Exceptions : throw on incorrect arguments
               throw on sql error (e.g. duplicate unique id)
  Caller     : general
  Status     : Stable

=cut
*/
/* NIY
sub store_alt_alleles {
  my $self  = shift;
  my $genes = shift;

  if (!ref($genes) eq 'ARRAY') {
    throw('List reference of Bio::EnsEMBL::Gene argument expected.');
  }

  my @genes     = @$genes;
  my $num_genes = scalar(@genes);

  if ($num_genes < 2) {
    warning('At least 2 genes must be provided to construct alternative alleles (gene id: ' . $genes[0]->dbID() . '). Ignoring.');
    return;
  }

  my @is_ref;
  my @ref_genes     = ();
  my @non_ref_genes = ();
  my @gene_ids      = ();

  foreach my $gene (@genes) {

    if (!ref($gene) || !$gene->isa('Bio::EnsEMBL::Gene')) {
      throw('List reference of Bio::EnsEMBL::Gene argument expected.');
    }

    my $gene_id = $gene->dbID();

    if (!$gene_id) {
      throw('Genes must have dbIDs in order to construct alternative alleles.');
    } else {
      push @gene_ids, $gene_id;
    }

    my $is_ref = $gene->slice->is_reference();

    push @is_ref, $is_ref;

    if ($is_ref) {
      push @ref_genes, $gene->dbID();
    } else {
      push @non_ref_genes, $gene->dbID();
    }
  }
  if (scalar(@ref_genes) > 1) {
    warning('More than one alternative allele on the reference sequence (gene ids: ' . join(',', @ref_genes) . '). Ignoring.');
    return;
  }

  #
  #insert the first gene seperately in order to get a unique identifier for
  #the set of alleles
  #

  my $sth = $self->prepare("INSERT INTO alt_allele (gene_id, is_ref) VALUES (?,?)");
  $sth->bind_param(1, $gene_ids[0], SQL_INTEGER);
  $sth->bind_param(2, $is_ref[0],   SQL_INTEGER);
  eval { $sth->execute(); };
  my $alt_allele_id = $sth->{'mysql_insertid'};

  if (!$alt_allele_id || $@) {
    throw("An SQL error occured inserting alternative alleles:\n$@");
  }
  $sth->finish();
  #
  # Insert all subsequent alt alleles using the alt_allele identifier
  # from the first insert
  #

  $sth = $self->prepare("INSERT INTO alt_allele (alt_allele_id, gene_id, is_ref) " . "VALUES (?,?,?)");

  for (my $i = 1; $i < $num_genes; $i++) {

    $sth->bind_param(1, $alt_allele_id, SQL_INTEGER);
    $sth->bind_param(2, $gene_ids[$i],  SQL_INTEGER);
    $sth->bind_param(3, $is_ref[$i],    SQL_INTEGER);
    eval { $sth->execute(); };

    if ($@) {
      # an error occured, revert the db to the previous state
      $sth = $self->prepare("DELETE FROM alt_allele WHERE alt_allele_id = ?");
      $sth->bind_param(1, $alt_allele_id, SQL_INTEGER);
      $sth->execute();
      $sth->finish();
      throw("An SQL error occured inserting alternative alleles:\n$@");
    }
  }

  $sth->finish();

  return $alt_allele_id;
} ## end sub store_alt_alleles
*/

/*
=head2 store

  Arg [1]    : Bio::EnsEMBL::Gene $gene
               The gene to store in the database
  Arg [2]    : ignore_release in xrefs [default 1] set to 0 to use release info 
               in external database references
  Example    : $gene_adaptor->store($gene);
  Description: Stores a gene in the database.
  Returntype : the database identifier (dbID) of the newly stored gene
  Exceptions : thrown if the $gene is not a Bio::EnsEMBL::Gene or if 
               $gene does not have an analysis object
  Caller     : general
  Status     : Stable

=cut
*/
// Default should be ignoreRelease = 1 so if you don't care about ignoreRelease
// use that
IDType GeneAdaptor_store(GeneAdaptor *ga, Gene *gene, int ignoreRelease)  {
  int i;

  if (gene == NULL) {
    fprintf(stderr, "feature is NULL in Gene_store\n");
    exit(1);
  }

  Class_assertType(CLASS_GENE, gene->objectType);

/* Should always be defined - its an arg to method
  if (!defined($ignore_release)) {
    $ignore_release = 1;
  }
*/
  DBAdaptor *db = ga->dba;
  AnalysisAdaptor *analysisAdaptor = DBAdaptor_getAnalysisAdaptor(db);

  if (Gene_isStored(gene, db)) {
    fprintf(stderr, "Gene ["IDFMTSTR"] is already stored in this database.\n", Gene_getDbID(gene) );
    return Gene_getDbID(gene);
  }


  // ensure coords are correct before storing
/* NIY!!!!!!!!!!!!!
  $gene->recalculate_coordinates();
*/

  Analysis *analysis = Gene_getAnalysis(gene);
  if (analysis == NULL) {
    fprintf(stderr,"An analysis must be attached to the features to be stored.\n");
    exit(1);
  }

  // store the analysis if it has not been stored yet
  if (!Analysis_isStored(analysis, db)) {
    fprintf(stderr, "STORING ANALYSIS %s\n", Analysis_getLogicName(analysis));
    AnalysisAdaptor_store(analysisAdaptor, analysis);
  }
  IDType analysisId = Analysis_getDbID(analysis);

/* Think above is equivalent to this
  if ($analysis->is_stored($db)) {
    $analysis_id = $analysis->dbID();
  } else {
    $analysis_id = $db->get_AnalysisAdaptor->store($analysis);
  }
*/

  char *type;
  if (Gene_getBiotype(gene)) {
    type = Gene_getBiotype(gene);
  } else {
    type = emptyString; /* static "" */
  }

  // default to is_current = 1 if this attribute is not set
// In C I initialise this to 1 in Gene_new
  int isCurrent = Gene_getIsCurrent(gene);
/*
  my $is_current = $gene->is_current;
  $is_current = 1 unless (defined($is_current));
*/

/* Note not doing transfer currently. If this is necessary?? then I'll have to rework this
  my $original             = $gene;
  my $original_transcripts = $gene->get_all_Transcripts();
  my $seq_region_id;

  ($gene, $seq_region_id) = $self->_pre_store($gene);
*/

  IDType seqRegionId = BaseFeatureAdaptor_preStore((BaseFeatureAdaptor *)ga, (SeqFeature*)gene); 

  char fmtStr[1024];
  char qStr[1024];
 
  // Canonical transcript ID will be updated later.
  // Set it to zero for now.
  sprintf(fmtStr, 
        "INSERT INTO gene "
           "SET analysis_id = "IDFMTSTR","
               "seq_region_id = "IDFMTSTR","
               "seq_region_start = %ld,"
               "seq_region_end = %ld,"
               "seq_region_strand = %d,"
               "biotype = '%s',"
               "description = %%s,"
               "source = %%s,"
               "status = %%s,"
               //"description = '%s',"
               //"source = '%s',"
               //"status = '%s',"
               "is_current = %d,"
               "canonical_transcript_id = 0," ,
               //"canonical_annotation = '%s'", 
               //"canonical_annotation = %%s", 
           analysisId, 
           seqRegionId, 
          Gene_getSeqRegionStart((SeqFeature*)gene),
          Gene_getSeqRegionEnd((SeqFeature*)gene),
          Gene_getSeqRegionStrand((SeqFeature*)gene),
           type, 
           //Gene_getDescription(gene),
           //Gene_getSource(gene),
           //Gene_getStatus(gene),
           isCurrent);
           //Gene_getCanonicalAnnotation(gene));

  char descQStr[1024];
  Gene_getDescription(gene) ? sprintf(descQStr,"'%s'", Gene_getDescription(gene)) : sprintf(descQStr, "NULL");

  char sourceQStr[1024];
  Gene_getSource(gene) ? sprintf(sourceQStr,"'%s'", Gene_getSource(gene)) : sprintf(sourceQStr, "NULL");

  char statusQStr[1024];
  Gene_getStatus(gene) ? sprintf(statusQStr,"'%s'", Gene_getStatus(gene)) : sprintf(statusQStr, "NULL");

  char canAnnQStr[1024];
  Gene_getCanonicalAnnotation(gene) ? sprintf(canAnnQStr,"'%s'", Gene_getCanonicalAnnotation(gene)) : sprintf(canAnnQStr, "NULL");

  sprintf(qStr, fmtStr, descQStr, sourceQStr, statusQStr, canAnnQStr);

  if (Gene_getStableId(gene)) {
/* Use FROM_UNIXTIME for now
    my $created  = $self->db->dbc->from_seconds_to_date($gene->created_date());
    my $modified = $self->db->dbc->from_seconds_to_date($gene->modified_date());
*/
  
    int version = Gene_getVersion(gene) > 0 ? Gene_getVersion(gene) : 1; // Assume version will be positive, Gene sets it to -1 when initialised
    sprintf(qStr,"%s, stable_id = '%s', version = %d, created_date = FROM_UNIXTIME(%ld), modified_date = FROM_UNIXTIME(%ld)",
            qStr, Gene_getStableId(gene), version, Gene_getCreated(gene), Gene_getModified(gene));
  }

  StatementHandle *sth = ga->prepare((BaseAdaptor *)ga,qStr,strlen(qStr));

  sth->execute(sth);
// NIY??? Finish was before insert id call?? Moved before
  IDType geneId = sth->getInsertId(sth);

  sth->finish(sth);


  // store the dbentries associated with this gene
/* NIY
  DBEntryAdaptor *dbEntryAdaptor = DBAdaptor_getDBEntryAdaptor(db);

  Vector *dbEntries = Gene_getAllDBEntries(gene);
  int i;
  for (i=0; i<Vector_getNumElement(dbEntries); i++) {
    DBEntry *dbe = Vector_getElementAt(dbEntries, i);
    DBEntryAdaptor_store(dbEntryAdaptor, dbe, geneId, "Gene", ignoreRelease);
  }
*/

  // We allow transcripts not to share equal exons and instead have
  // copies.  For the database we still want sharing though, to have
  // easier time with stable ids. So we need to have a step to merge
  // exons together before store.
  StringHash *exonHash = StringHash_new(STRINGHASH_SMALL);

  for (i=0; i<Gene_getTranscriptCount(gene); i++) {
    Transcript *trans = Gene_getTranscriptAt(gene, i);
    int j;
    for (j=0; j<Transcript_getExonCount(trans); j++) {
      Exon *exon = Transcript_getExonAt(trans, j);
      char exonKey[2048];
      Exon_getHashKey(exon, exonKey);

      if (StringHash_contains(exonHash, exonKey)) {
        Transcript_swapExons(trans, exon, (Exon *)StringHash_getValue(exonHash, exonKey));
      } else {
        StringHash_add(exonHash, exonKey, exon);
      }
    }
  }
  StringHash_free(exonHash, NULL);

  TranscriptAdaptor *transcriptAdaptor = DBAdaptor_getTranscriptAdaptor(db);

//  Vector *transcripts = Gene_getAllTranscripts(gene);
//  my $transcripts = $gene->get_all_Transcripts();

  IDType newCanonicalTranscriptId = 0; // Assume 0 is not a valid dbID - NIY May not be a good assumption

// Note I'm not doing the transfer so no new/old separation for now.
// Just call it new
  for (i = 0; i < Gene_getTranscriptCount(gene); i++) {
    Transcript *new = Gene_getTranscriptAt(gene, i);
//    my $old = $original_transcripts->[$i];

    TranscriptAdaptor_store(transcriptAdaptor, new, geneId, analysisId);

    if (!newCanonicalTranscriptId && Transcript_getIsCanonical(new)) {
      newCanonicalTranscriptId = Transcript_getDbID(new);
    }

    // update the original transcripts since we may have made copies of
    // them by transforming the gene
/*
    $old->dbID($new->dbID());
    $old->adaptor($new->adaptor());

    if ($new->translation) {
      $old->translation->dbID($new->translation()->dbID);
      $old->translation->adaptor($new->translation()->adaptor);
    }
*/
  }

//  if (defined($new_canonical_transcript_id)) {
  if (newCanonicalTranscriptId) {
    // Now the canonical transcript has been stored, so update the
    // canonical_transcript_id of this gene with the new dbID.
    sprintf(qStr,"UPDATE gene SET canonical_transcript_id = "IDFMTSTR" WHERE gene_id = "IDFMTSTR,newCanonicalTranscriptId, geneId);

    sth = ga->prepare((BaseAdaptor *)ga,qStr,strlen(qStr));

    sth->execute(sth);
    sth->finish(sth);
  }

  // update gene to point to display xref if it is set
/* NIY
  DBEntry *displayXref = Gene_getDisplayXref(gene);
  if (displayXref != NULL) {

    IDType dxrefId = 0;
    if (DBEntry_isStored(displayXref, db)) {
      dxrefId = DBEntry_getDbID(displayXref);
    } else {
      dxrefId = DBEntryAdaptor_exists(dbEntryAdaptor, displayXref);
    }

//    if (defined($dxref_id)) {
    if (dxrefId) {
      sprintf(qStr, "UPDATE gene SET display_xref_id = "IDFMTSTR" WHERE gene_id = "IDFMTSTR, dxrefId, geneId);

      StatementHandle *sth = ga->prepare((BaseAdaptor *)ga,qStr,strlen(qStr));

      sth->execute(sth);
      sth->finish(sth);

      DBEntry_setDbID(displayXref, dxrefId);
      DBEntry_setAdaptor(displayXref, (BaseAdaptor *)dbEntryAdaptor);
// Duplicate code in perl
//      $display_xref->dbID($dxref_id);
//      $display_xref->adaptor($dbEntryAdaptor);
    } else {
      fprintf(stderr, "Display_xref %s:%s is not stored in database.\n"
                      "Not storing relationship to this gene.\n", 
              DBEntry_getDbName(displayXref), DBEntry_getDisplayId(displayXref));
      DBEntry_setDbID(displayXref, 0);
      DBEntry_setAdaptor(displayXref, NULL);
    }
  }
*/

  // store gene attributes if there are any
  AttributeAdaptor *attrAdaptor = DBAdaptor_getAttributeAdaptor(db);
  Vector *attribs = Gene_getAllAttributes(gene, NULL);
  AttributeAdaptor_storeOnGeneId(attrAdaptor, geneId, attribs);
  Vector_free(attribs);

  // store unconventional transcript associations if there are any
/* NIY - and probably never will be implemented, not used anywhere that I know of
  my $utaa = $db->get_UnconventionalTranscriptAssociationAdaptor();
  foreach my $uta (@{$gene->get_all_unconventional_transcript_associations()}) {
    $utaa->store($uta);
  }
*/

  // set the adaptor and dbID on the original passed in gene not the
  // transfered copy - I don't make the transferred copy in C (hopefully I don't need to)
  Gene_setAdaptor(gene, (BaseAdaptor *)ga);
  Gene_setDbID(gene,geneId);

  return Gene_getDbID(gene);
}

/* Don't bother
=head2 remove

  Arg [1]    : Bio::EnsEMBL::Gene $gene
               the gene to remove from the database
  Example    : $gene_adaptor->remove($gene);
  Description: Removes a gene completely from the database. All associated
               transcripts, exons, stable_identifiers, descriptions, etc.
               are removed as well. Use with caution!
  Returntype : none
  Exceptions : throw on incorrect arguments 
               warning if gene is not stored in this database
  Caller     : general
  Status     : Stable

=cut

sub remove {
  my $self = shift;
  my $gene = shift;

  if (!ref($gene) || !$gene->isa('Bio::EnsEMBL::Gene')) {
    throw("Bio::EnsEMBL::Gene argument expected.");
  }

  if (!$gene->is_stored($self->db())) {
    warning("Cannot remove gene " . $gene->dbID() . ". Is not stored in " . "this database.");
    return;
  }

  # remove all object xrefs associated with this gene

  my $dbe_adaptor = $self->db()->get_DBEntryAdaptor();
  foreach my $dbe (@{$gene->get_all_DBEntries()}) {
    $dbe_adaptor->remove_from_object($dbe, $gene, 'Gene');
  }

  # remove all alternative allele entries associated with this gene
  my $sth = $self->prepare("DELETE FROM alt_allele WHERE gene_id = ?");
  $sth->bind_param(1, $gene->dbID, SQL_INTEGER);
  $sth->execute();
  $sth->finish();

  # remove the attributes associated with this transcript
  my $attrib_adaptor = $self->db->get_AttributeAdaptor;
  $attrib_adaptor->remove_from_Gene($gene);

  # remove all of the transcripts associated with this gene
  my $transcriptAdaptor = $self->db->get_TranscriptAdaptor();
  foreach my $trans (@{$gene->get_all_Transcripts()}) {
    $transcriptAdaptor->remove($trans);
  }

  # remove any unconventional transcript associations involving this gene

  $sth = $self->prepare("DELETE FROM unconventional_transcript_association " . "WHERE gene_id = ? ");
  $sth->bind_param(1, $gene->dbID, SQL_INTEGER);
  $sth->execute();
  $sth->finish();

  # remove this gene from the database

  $sth = $self->prepare("DELETE FROM gene WHERE gene_id = ? ");
  $sth->bind_param(1, $gene->dbID, SQL_INTEGER);
  $sth->execute();
  $sth->finish();

  # unset the gene identifier and adaptor thereby flagging it as unstored

  $gene->dbID(undef);
  $gene->adaptor(undef);

  return;
} ## end sub remove
*/

/* Don't bother
=head2 get_Interpro_by_geneid

  Arg [1]    : String $gene_stable_id
               The stable ID of the gene to obtain
  Example    : @i = @{
                  $gene_adaptor->get_Interpro_by_geneid(
                    $gene->stable_id() ) };
  Description: Gets interpro accession numbers by gene stable id. A hack really
               - we should have a much more structured system than this.
  Returntype : listref of strings (Interpro_acc:description)
  Exceptions : none 
  Caller     : domainview
  Status     : Stable

=cut

sub get_Interpro_by_geneid {
  my ($self, $gene_stable_id) = @_;

  my $sql = qq(
  SELECT    i.interpro_ac,
            x.description
  FROM      transcript t,
            translation tl,
            protein_feature pf,
            interpro i,
            xref x,
            gene g
  WHERE     g.stable_id = ?
    AND     t.gene_id = g.gene_id
    AND     t.is_current = 1
    AND     tl.transcript_id = t.transcript_id
    AND     tl.translation_id = pf.translation_id
    AND     i.id = pf.hit_name
    AND     i.interpro_ac = x.dbprimary_acc);

  my $sth = $self->prepare($sql);

  $sth->bind_param(1, $gene_stable_id, SQL_VARCHAR);

  $sth->execute;

  my @out;
  my %h;
  while ((my $arr = $sth->fetchrow_arrayref())) {
    if ($h{$arr->[0]}) { next; }
    $h{$arr->[0]} = 1;
    my $string = $arr->[0] . ":" . $arr->[1];
    push(@out, $string);
  }

  return \@out;
} ## end sub get_Interpro_by_geneid
*/

/*
=head2 update

  Arg [1]    : Bio::EnsEMBL::Gene $gene
               The gene to update
  Example    : $gene_adaptor->update($gene);
  Description: Updates the type, analysis, display_xref, status, is_current and
               description of a gene in the database.
  Returntype : None
  Exceptions : thrown if the $gene is not a Bio::EnsEMBL::Gene
  Caller     : general
  Status     : Stable

=cut
*/
/* NIY
sub update {
  my ($self, $gene) = @_;
  my $update = 0;

  if (!defined $gene || !ref $gene || !$gene->isa('Bio::EnsEMBL::Gene')) {
    throw("Must update a gene object, not a $gene");
  }

  my $update_gene_sql = qq(
       UPDATE gene
          SET biotype = ?,
              analysis_id = ?,
              display_xref_id = ?,
              status = ?,
              description = ?,
              is_current = ?,
              canonical_transcript_id = ?,
              canonical_annotation = ?
        WHERE gene_id = ?
  );

  my $display_xref = $gene->display_xref();
  my $display_xref_id;

  if ($display_xref && $display_xref->dbID()) {
    $display_xref_id = $display_xref->dbID();
  } else {
    $display_xref_id = undef;
  }

  my $sth = $self->prepare($update_gene_sql);

  $sth->bind_param(1, $gene->biotype(),        SQL_VARCHAR);
  $sth->bind_param(2, $gene->analysis->dbID(), SQL_INTEGER);
  $sth->bind_param(3, $display_xref_id,        SQL_INTEGER);
  $sth->bind_param(4, $gene->status(),         SQL_VARCHAR);
  $sth->bind_param(5, $gene->description(),    SQL_VARCHAR);
  $sth->bind_param(6, $gene->is_current(),     SQL_TINYINT);

  if (defined($gene->canonical_transcript())) {
    $sth->bind_param(7, $gene->canonical_transcript()->dbID(), SQL_INTEGER);
  } else {
    $sth->bind_param(7, 0, SQL_INTEGER);
  }

  $sth->bind_param(8, $gene->canonical_annotation(), SQL_VARCHAR);
  $sth->bind_param(9, $gene->dbID(), SQL_INTEGER);

  $sth->execute();

  # maybe should update stable id ???
} ## end sub update
*/

/*
# _objs_from_sth

#  Arg [1]    : StatementHandle $sth
#  Arg [2]    : Bio::EnsEMBL::AssemblyMapper $mapper
#  Arg [3]    : Bio::EnsEMBL::Slice $dest_slice
#  Description: PROTECTED implementation of abstract superclass method.
#               responsible for the creation of Genes
#  Returntype : listref of Bio::EnsEMBL::Genes in target coordinate system
#  Exceptions : none
#  Caller     : internal
#  Status     : Stable
*/
Vector *GeneAdaptor_objectsFromStatementHandle(GeneAdaptor *ga,
                                               StatementHandle *sth,
                                               AssemblyMapper *assMapper,
                                               Slice *destSlice) {
  // This code is ugly because an attempt has been made to remove as many
  // function calls as possible for speed purposes.  Thus many caches and
  // a fair bit of gymnastics is used.
  //
  SliceAdaptor *sa     = DBAdaptor_getSliceAdaptor(ga->dba);
  AnalysisAdaptor *aa  = DBAdaptor_getAnalysisAdaptor(ga->dba);
  DBEntryAdaptor *dbea = DBAdaptor_getDBEntryAdaptor(ga->dba);

  Vector *genes = Vector_new();
  IDHash *sliceHash = IDHash_new(IDHASH_SMALL);
/* Don't bother with these three - analysis is cached in its adaptor, and I can't believe speed 
   is going to be limited by name and cs access functions on a slice!
  my %analysis_hash;
  my %sr_name_hash;
  my %sr_cs_hash;
*/


/* Unused!
  CoordSystem *asmCs;
  CoordSystem *cmpCs;
  my $asm_cs_vers;
  my $asm_cs_name;
  my $cmp_cs_vers;
  my $cmp_cs_name;

  if ($mapper) {
    $asm_cs      = $mapper->assembled_CoordSystem();
    $cmp_cs      = $mapper->component_CoordSystem();
    $asm_cs_name = $asm_cs->name();
    $asm_cs_vers = $asm_cs->version();
    $cmp_cs_name = $cmp_cs->name();
    $cmp_cs_vers = $cmp_cs->version();
  }
*/

  long   destSliceStart;
  long   destSliceEnd;
  int    destSliceStrand;
  long   destSliceLength;
  char * destSliceSrName;
  IDType destSliceSrId = 0;

  if (destSlice) {
    destSliceStart  = Slice_getStart(destSlice);
    destSliceEnd    = Slice_getEnd(destSlice);
    destSliceStrand = Slice_getStrand(destSlice);
    destSliceLength = Slice_getLength(destSlice);
    destSliceSrName = Slice_getSeqRegionName(destSlice);
    destSliceSrId   = Slice_getSeqRegionId(destSlice);
  }

// Note FEATURE label is here
//FEATURE: while ($sth->fetch()) 
  ResultRow *row;
  while ((row = sth->fetchRow(sth))) {
    IDType geneId =                row->getLongLongAt(row,0);
    IDType seqRegionId =           row->getLongLongAt(row,1);
    long seqRegionStart =          row->getLongAt(row,2);
    long seqRegionEnd =            row->getLongAt(row,3);
    int seqRegionStrand =          row->getIntAt(row,4);
    IDType analysisId =            row->getLongLongAt(row,5);
    char *biotype =                row->getStringAt(row,6);
    IDType displayXrefId =         row->getLongLongAt(row,7);
    char *geneDescription =        row->getStringAt(row,8);
    char *status =                 row->getStringAt(row,9);
    char *source =                 row->getStringAt(row,10);
    int isCurrent =                row->getIntAt(row,11);
    IDType canonicalTranscriptId = row->getLongLongAt(row,12);
    char *canonicalAnnotation =    row->getStringAt(row,13);
    char *stableId =               row->getStringAt(row,14);
    int version =                  row->getIntAt(row,15);
    int createdDate =              row->getIntAt(row,16);
    int modifiedDate =             row->getIntAt(row,17);
// Changed from perl xrefDisplayId  to xrefDisplayLabel to match Transcript version
    char *xrefDisplayLabel =       row->getStringAt(row,18);
    char *xrefPrimaryAcc =         row->getStringAt(row,19);
    char *xrefDesc =               row->getStringAt(row,20);
    char *xrefVersion =            row->getStringAt(row,21);
    char *externalDb =             row->getStringAt(row,22);
    char *externalStatus =         row->getStringAt(row,23);
    char *externalRelease =        row->getStringAt(row,24);
    char *externalDbName =         row->getStringAt(row,25);
// Added xref prefix to info* var names to match Transcript version
    char *xrefInfoType =           row->getStringAt(row,26);
    char *xrefInfoText =           row->getStringAt(row,27);

    // get the analysis object
    Analysis *analysis = AnalysisAdaptor_fetchByDbID(aa, analysisId);

    // Not doing internal seq id stuff for now
    // need to get the internal_seq_region, if present
    //$seq_region_id = $self->get_seq_region_id_internal($seq_region_id);
    //my $slice = $slice_hash{"ID:" . $seq_region_id};
    //if (!$slice) {
    //  $slice                              = $sa->fetch_by_seq_region_id($seq_region_id);
    //  $slice_hash{"ID:" . $seq_region_id} = $slice;
    // // $sr_name_hash{$seq_region_id}       = $slice->seq_region_name();
    // //  $sr_cs_hash{$seq_region_id}         = $slice->coord_system();
    //}

    if (! IDHash_contains(sliceHash, seqRegionId)) {
      IDHash_add(sliceHash, seqRegionId, SliceAdaptor_fetchBySeqRegionId(sa, seqRegionId, POS_UNDEF, POS_UNDEF, STRAND_UNDEF));
    }
    Slice *slice = IDHash_getValue(sliceHash, seqRegionId);

    Slice *geneSlice = slice;

    char *srName      = Slice_getSeqRegionName(slice);
    CoordSystem *srCs = Slice_getCoordSystem(slice);

    //
    // remap the feature coordinates to another coord system
    // if a mapper was provided
    //
    if (assMapper != NULL) {
      MapperRangeSet *mrs;

      // Slightly suspicious about need for this if statement so left in perl statements for now
      if (destSlice != NULL &&
          assMapper->objectType == CLASS_CHAINEDASSEMBLYMAPPER) {
        mrs = ChainedAssemblyMapper_map(assMapper, srName, seqRegionStart, seqRegionEnd, seqRegionStrand, srCs, 1, destSlice);

        //($seq_region_id, $seq_region_start, $seq_region_end, $seq_region_strand) = $mapper->map($sr_name, $seq_region_start, $seq_region_end, $seq_region_strand, $sr_cs, 1, $dest_slice);

      } else {
        mrs = AssemblyMapper_fastMap(assMapper, srName, seqRegionStart, seqRegionEnd, seqRegionStrand, srCs, NULL);
        //($seq_region_id, $seq_region_start, $seq_region_end, $seq_region_strand) = $mapper->fastmap($sr_name, $seq_region_start, $seq_region_end, $seq_region_strand, $sr_cs);
      }

      // skip features that map to gaps or coord system boundaries
      //next FEATURE if (!defined($seq_region_id));
      if (MapperRangeSet_getNumRange(mrs) == 0) {
        continue;
      }
      MapperRange *range = MapperRangeSet_getRangeAt(mrs, 0);
      if (range->rangeType == MAPPERRANGE_GAP) {
        fprintf(stderr,"Got a mapper gap in gene obj_from_sth - not sure if this is allowed\n");
        exit(1);
      } else {
        MapperCoordinate *mc = (MapperCoordinate *)range;

        seqRegionId     = mc->id;
        seqRegionStart  = mc->start;
        seqRegionEnd    = mc->end;
        seqRegionStrand = mc->strand;
      }

      MapperRangeSet_free(mrs);

      // get a slice in the coord system we just mapped to
      //$slice = $slice_hash{"ID:" . $seq_region_id} ||= $sa->fetch_by_seq_region_id($seq_region_id);
      if (! IDHash_contains(sliceHash, seqRegionId)) {
        IDHash_add(sliceHash, seqRegionId, SliceAdaptor_fetchBySeqRegionId(sa, seqRegionId, POS_UNDEF, POS_UNDEF, STRAND_UNDEF));
      }
      geneSlice = IDHash_getValue(sliceHash, seqRegionId);
    }

    // 
    // If a destination slice was provided convert the coords.
    //
    if (destSlice != NULL) {
      if (destSliceStrand == 1) {
        // Positive strand.
        seqRegionStart = seqRegionStart - destSliceStart + 1;
        seqRegionEnd   = seqRegionEnd - destSliceStart + 1;

// No circular stuff
        if (0) {
//        if ($dest_slice->is_circular())
/*
          # Handle cicular chromosomes.

          if ($seq_region_start > $seq_region_end) {
            # Looking at a feature overlapping the chromsome origin.

            if ($seq_region_end > $dest_slice_start) {
              # Looking at the region in the beginning of the
              # chromosome.
              $seq_region_start -= $dest_slice->seq_region_length();
            }

            if ($seq_region_end < 0) {
              $seq_region_end += $dest_slice->seq_region_length();
            }
          } else {
            if (destSliceStart > destSliceEnd && seqRegionEnd < 0) {
              // Looking at the region overlapping the chromosome
              // origin and a feature which is at the beginning of the
              // chromosome.
// Perl actually does two length calls here - I'd have thought this would be more expensive than the coordsys accessing call above!
              seqRegionStart += Slice_getSeqRegionLength(destSlice);
              seqRegionEnd   += Slice_getSeqRegionLength(destSlice);
            }
          }
*/
        }
      } else {
        // Negative strand.
// No circular stuff
        if (0) {
//        if (   $dest_slice->is_circular()
//            && $seq_region_start > $seq_region_end)
/*
          # Handle cicular chromosomes.

          if ($seq_region_end > $dest_slice_start) {
            # Looking at the region in the beginning of the
            # chromosome.
            $seq_region_start = $dest_slice_end - $seq_region_end + 1;
            $seq_region_end   = $seq_region_end - $dest_slice->seq_region_length - $dest_slice_start + 1;
          } else {
            my $tmp_seq_region_start = $seq_region_start;
            $seq_region_start = $dest_slice_end - $seq_region_end - $dest_slice->seq_region_length + 1;
            $seq_region_end   = $dest_slice_end - $tmp_seq_region_start + 1;
          }

*/
        } else {
          // Non-circular chromosome - sanity
          long tmpSeqRegionStart = seqRegionStart;
          seqRegionStart = destSliceEnd - seqRegionEnd + 1;
          seqRegionEnd   = destSliceEnd - tmpSeqRegionStart + 1;
        }

        seqRegionStrand = -seqRegionStrand;
      }

      // Throw away features off the end of the requested slice or on
      // different seq_region.
// Perl used 'ne' for comparison of ids but these are ints so that's not really very efficient
      if (seqRegionEnd < 1 || 
          seqRegionStart > destSliceLength || 
          (destSliceSrId != seqRegionId)) {
// Any freeing to do - don't think so??
        //next FEATURE;
        continue;
      }
      geneSlice = destSlice;
    }

    DBEntry *displayXref = NULL;

    // I presume displayXrefId will be zero if left join doesn't find one
    if (displayXrefId) {
      displayXref = DBEntry_new();

// Slightly wierd formatting to make it easier to read
      DBEntry_setAdaptor    (displayXref,   (BaseAdaptor *)dbea);
      DBEntry_setDbID       (displayXref,   displayXrefId);
      DBEntry_setPrimaryId  (displayXref,   xrefPrimaryAcc);
      DBEntry_setDisplayId  (displayXref,   xrefDisplayLabel);
      DBEntry_setVersion    (displayXref,   xrefVersion);
      DBEntry_setDescription(displayXref,   xrefDesc);
      DBEntry_setRelease    (displayXref,   externalRelease);
      DBEntry_setDbName     (displayXref,   externalDb);
      DBEntry_setDbDisplayName(displayXref, externalDbName);
      DBEntry_setInfoType   (displayXref,   xrefInfoType);
      DBEntry_setInfoText   (displayXref,   xrefInfoText);
      DBEntry_setStatus     (displayXref,   externalStatus);
    }

    // Finally, create the new Gene.
    Gene *gene = Gene_new();
    Gene_setAnalysis      (gene, analysis);
    Gene_setBiotype       (gene, biotype);
    Gene_setStart         (gene, seqRegionStart);
    Gene_setEnd           (gene, seqRegionEnd);
    Gene_setStrand        (gene, seqRegionStrand);
    Gene_setAdaptor       (gene, (BaseAdaptor *)ga);
    Gene_setSlice         (gene, geneSlice);
    Gene_setDbID          (gene, geneId);
    Gene_setStableId      (gene, stableId);
    Gene_setVersion       (gene, version);
// Had $created_date || undef, for this (and equivalent for modified_date - not sure about that???
    Gene_setCreated       (gene, createdDate);
    Gene_setModified      (gene, modifiedDate);
    Gene_setDescription   (gene, geneDescription);
    Gene_setExternalName  (gene, NULL);
    Gene_setExternalDb    (gene, externalDb);
    Gene_setExternalStatus(gene, externalStatus);
    Gene_setDisplayXref   (gene, displayXref);
    Gene_setStatus        (gene, status);
    Gene_setSource        (gene, source);
    Gene_setIsCurrent     (gene, isCurrent);
    Gene_setCanonicalTranscriptId(gene, canonicalTranscriptId);
    Gene_setCanonicalAnnotation(gene, canonicalAnnotation);

    Vector_addElement(genes, gene);
  }

  // Don't free slices because they might be being used????
  IDHash_free(sliceHash, NULL);

  return genes;
}

/*
=head2 cache_gene_seq_mappings

  Example    : $gene_adaptor->cache_gene_seq_mappings();
  Description: caches all the assembly mappings needed for genes
  Returntype : None
  Exceptions : None
  Caller     : general
  Status     : At Risk
             : New experimental code

=cut
*/
void GeneAdaptor_cacheGeneSeqMappings(GeneAdaptor *ga) {
  // get the sequence level to map too
  char qStr[1024];
  sprintf(qStr, "SELECT name FROM coord_system WHERE attrib like '%%%%sequence_level%%%%' AND species_id = "IDFMTSTR, 
          GeneAdaptor_getSpeciesId(ga));

  StatementHandle *sth = ga->prepare((BaseAdaptor *)ga,qStr,strlen(qStr));
  sth->execute(sth);

  ResultRow *row = sth->fetchRow(sth);
  char *sequenceLevel = row->getStringCopyAt(row, 0);

  sth->finish(sth);

  CoordSystemAdaptor *csa    = DBAdaptor_getCoordSystemAdaptor(ga->dba);
  AssemblyMapperAdaptor *ama = DBAdaptor_getAssemblyMapperAdaptor(ga->dba);

  CoordSystem *cs1 = CoordSystemAdaptor_fetchByName(csa, sequenceLevel, NULL);

  // get levels to map to 
  MetaCoordContainer *mcc   = DBAdaptor_getMetaCoordContainer(ga->dba);
  Vector *csNew = MetaCoordContainer_fetchAllCoordSystemsByFeatureType(mcc, "gene");

  int i;
  for (i=0; i<Vector_getNumElement(csNew); i++) {
    CoordSystem *cs2 = Vector_getElementAt(csNew, i);

    AssemblyMapper *am = AssemblyMapperAdaptor_fetchByCoordSystems(ama, cs1, cs2);
    AssemblyMapper_registerAll(am);
  }
 
  free(sequenceLevel);

  Vector_free(csNew);
  return;
}

/* Don't bother
=head2 fetch_all_by_exon_supporting_evidence

  Arg [1]    : String $hit_name
               Name of supporting feature
  Arg [2]    : String $feature_type 
               one of "dna_align_feature" or "protein_align_feature"
  Arg [3]    : (optional) Bio::Ensembl::Analysis
  Example    : $genes = $gene_adaptor->fetch_all_by_exon_supporting_evidence(
                  'XYZ', 'dna_align_feature');
  Description: Gets all the genes with transcripts with exons which have a
               specified hit on a particular type of feature. Optionally filter
               by analysis.
  Returntype : Listref of Bio::EnsEMBL::Gene
  Exceptions : If feature_type is not of correct type.
  Caller     : general
  Status     : Stable

=cut
sub fetch_all_by_exon_supporting_evidence {
  my ($self, $hit_name, $feature_type, $analysis) = @_;

  if ($feature_type !~ /(dna)|(protein)_align_feature/) {
    throw("feature type must be dna_align_feature or protein_align_feature");
  }

  my $anal_from = ", analysis a " if ($analysis);
  my $anal_where = "AND a.analysis_id = f.analysis_id AND a.analysis_id=? "
    if ($analysis);

  my $sql = qq(
      SELECT DISTINCT(g.gene_id)
        FROM gene g,
             transcript t,
             exon_transcript et,
             supporting_feature sf,
             $feature_type f
             $anal_from
       WHERE g.gene_id = t.gene_id
         AND g.is_current = 1
         AND t.transcript_id = et.transcript_id
         AND et.exon_id = sf.exon_id
         AND sf.feature_id = f.${feature_type}_id
         AND sf.feature_type = ?
         AND f.hit_name=?
         $anal_where
  );

  my $sth = $self->prepare($sql);

  $sth->bind_param(1, $feature_type, SQL_VARCHAR);
  $sth->bind_param(2, $hit_name,     SQL_VARCHAR);
  $sth->bind_param(3, $analysis->dbID(), SQL_INTEGER) if ($analysis);

  $sth->execute();

  my @genes;

  while (my $id = $sth->fetchrow_array) {
    my $gene = $self->fetch_by_dbID($id);
    push(@genes, $gene) if $gene;
  }

  return \@genes;
} ## end sub fetch_all_by_exon_supporting_evidence
*/

/* Don't bother
=head2 fetch_all_by_transcript_supporting_evidence

  Arg [1]    : String $hit_name
               Name of supporting feature
  Arg [2]    : String $feature_type 
               one of "dna_align_feature" or "protein_align_feature"
  Arg [3]    : (optional) Bio::Ensembl::Analysis
  Example    : $genes = $gene_adaptor->fetch_all_by_transcript_supporting_evidence('XYZ', 'dna_align_feature');
  Description: Gets all the genes with transcripts with evidence for a
               specified hit on a particular type of feature. Optionally filter
               by analysis.
  Returntype : Listref of Bio::EnsEMBL::Gene.
  Exceptions : If feature_type is not of correct type.
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_transcript_supporting_evidence {
  my ($self, $hit_name, $feature_type, $analysis) = @_;

  if ($feature_type !~ /(dna)|(protein)_align_feature/) {
    throw("feature type must be dna_align_feature or protein_align_feature");
  }

  my $anal_from = ", analysis a " if ($analysis);
  my $anal_where = "AND a.analysis_id = f.analysis_id AND a.analysis_id=? "
    if ($analysis);

  my $sql = qq(
      SELECT DISTINCT(g.gene_id)
        FROM gene g,
             transcript t,
             transcript_supporting_feature sf,
             $feature_type f
             $anal_from
       WHERE g.gene_id = t.gene_id
         AND g.is_current = 1
         AND t.transcript_id = sf.transcript_id
         AND sf.feature_id = f.${feature_type}_id
         AND sf.feature_type = ?
         AND f.hit_name=?
         $anal_where
  );

  my $sth = $self->prepare($sql);

  $sth->bind_param(1, $feature_type, SQL_VARCHAR);
  $sth->bind_param(2, $hit_name,     SQL_VARCHAR);
  $sth->bind_param(3, $analysis->dbID(), SQL_INTEGER) if ($analysis);

  $sth->execute();

  my @genes;

  while (my $id = $sth->fetchrow_array) {
    my $gene = $self->fetch_by_dbID($id);
    push(@genes, $gene) if $gene;
  }

  return \@genes;
} ## end sub fetch_all_by_transcript_supporting_evidence
*/

/* Don't bother - why does this exist!!!!
=head2 fetch_nearest_Gene_by_Feature

  Arg [1]    : Feature object
  Example    : $genes = $gene_adaptor->fetch_nearest_Gene_by_Feature($feat);
  Description: Gets the nearest gene to the feature 
  Returntype : Listref of Bio::EnsEMBL::Gene, EMPTY list if no nearest
  Caller     : general
  Status     : UnStable

=cut

sub fetch_nearest_Gene_by_Feature {
  my $self = shift;
  my $feat = shift;

  my $stranded = shift;
  my $stream   = shift;    # 1 up stream -1 downstream
  my @genes;

  my $strand = $feat->strand;
  if (defined($stream) and !$strand) {
    warn("stream specified but feature has no strand so +ve strand will be used");
    $strand = 1;
  }
  my $min_dist = 999;
  my $gene_id  = 0;

  my $overlapping = $feat->get_overlapping_Genes();

  return $overlapping if (defined(@{$overlapping}[0]));

  my $seq_region_id = $feat->slice->adaptor->get_seq_region_id($feat->slice);
  my $start         = ($feat->start + $feat->slice->start) - 1;
  my $end           = ($feat->end + $feat->slice->start) - 1;

  my @gene_ids;
  if (!defined($stream) or $stream == 0) {

    my $sql1 = "select g.gene_id, (? - g.seq_region_end)  as 'dist' from gene g where ";
    if ($stranded) {
      $sql1 .= "g.seq_region_strand = " . $strand . " and ";
    }
    $sql1 .= "seq_region_id = ? and g.seq_region_end < ? order by dist limit 10";

    #
    # MAYBE set the result of prepare to be static in case lots of calls.
    #
    my $sql1_sth = $self->prepare($sql1) || die "Could not prepare $sql1";
    $sql1_sth->execute($start, $seq_region_id, $start)
      || die "Could not execute sql";
    $sql1_sth->bind_columns(\$gene_id, \$min_dist)
      || die "Could mot bin columns";

    my $last_dist = 99999999999999999;
    while ($sql1_sth->fetch()) {
      if ($min_dist <= $last_dist) {
        push @gene_ids, $gene_id;
        $last_dist = $min_dist;
      }
    }
    $sql1_sth->finish();

    my $sql2 = "select g.gene_id, (g.seq_region_start - ?)  as 'dist' from gene g  where ";
    if ($stranded) {
      $sql2 .= "g.seq_region_strand = " . $feat->strand . " and ";
    }
    $sql2 .= "seq_region_id = ? and g.seq_region_start > ? order by dist limit 10";

    my $sql2_sth = $self->prepare($sql2) || die "could not prepare $sql2";

    my ($tmp_min_dist, $tmp_gene_id);
    $sql2_sth->execute($end, $seq_region_id, $end)
      || die "Could not execute sql";
    $sql2_sth->bind_columns(\$tmp_gene_id, \$tmp_min_dist)
      || die "Could mot bin columns";
    my $first = 1;
    while ($sql2_sth->fetch()) {
      if ($tmp_min_dist <= $last_dist) {
        if ($first) {
          $first = 0;
          if ($tmp_min_dist < $last_dist) {
            @gene_ids = ();    #reset
          }
        }
        push @gene_ids, $tmp_gene_id;
        $last_dist = $tmp_min_dist;
      }
    }
    $sql2_sth->finish();

  } elsif (($stream*$strand) == 1) {
    my $sql1 = "select g.gene_id, (? - g.seq_region_end)  as 'dist' from gene g where ";
    if ($stranded) {
      $sql1 .= "g.seq_region_strand = " . $strand . " and ";
    }
    $sql1 .= "seq_region_id = ? and g.seq_region_end < ? order by dist limit 10";

    #
    # MAYBE set the result of prepare to be static in case lots of calls.
    #
    my $sql1_sth = $self->prepare($sql1) || die "Could not prepare $sql1";
    $sql1_sth->execute($start, $seq_region_id, $start)
      || die "Could not execute sql";
    $sql1_sth->bind_columns(\$gene_id, \$min_dist)
      || die "Could mot bin columns";

    my $last_dist;
    my $first = 1;
    while ($sql1_sth->fetch()) {
      if ($first) {
        $first = 0;
      } else {
        next if ($min_dist > $last_dist);
      }
      push @gene_ids, $gene_id;
      $last_dist = $min_dist;
    }
    $sql1_sth->finish();
  } elsif (($stream*$strand) == -1) {

    my $sql2 = "select g.gene_id, (g.seq_region_start - ?)  as 'dist' from gene g  where ";
    if ($stranded) {
      $sql2 .= "g.seq_region_strand = " . $feat->strand . " and ";
    }
    $sql2 .= "seq_region_id = ? and g.seq_region_start > ? order by dist limit 10";

    my $sql2_sth = $self->prepare($sql2) || die "could not prepare $sql2";

    my ($tmp_min_dist, $tmp_gene_id);
    $sql2_sth->execute($end, $seq_region_id, $end)
      || die "Could not execute sql";
    $sql2_sth->bind_columns(\$tmp_gene_id, \$tmp_min_dist)
      || die "Could mot bin columns";
    my $first = 1;
    my $last_dist;
    while ($sql2_sth->fetch()) {
      if ($first) {
        $first = 0;
      } else {
        next if ($tmp_min_dist > $last_dist);
      }
      push @gene_ids, $tmp_gene_id;
      $last_dist = $tmp_min_dist;
    }
    $sql2_sth->finish();
  } else {
    die "Invalid stream or strand must be -1, 0 or 1\n";
  }

  foreach my $gene_id (@gene_ids) {
    push @genes, $self->fetch_by_dbID($gene_id);
  }
  return \@genes;

} ## end sub fetch_nearest_Gene_by_Feature
*/

/* Hopefully won't need
##########################
#                        #
#  DEPRECATED METHODS    #
#                        #
##########################

=head2 fetch_by_maximum_DBLink

 DEPRECATED - use fetch_all_by_external_name instead

=cut

sub fetch_by_maximum_DBLink {
  my ($self, $external_id) = @_;

  deprecate("use fetch_all_by_external_name instead");

  my $genes = $self->fetch_all_by_external_name($external_id);

  my $biggest;
  my $max  = 0;
  my $size = scalar(@$genes);
  if ($size > 0) {
    foreach my $gene (@$genes) {
      my $size = scalar(@{$gene->get_all_Exons});
      if ($size > $max) {
        $biggest = $gene;
        $max     = $size;
      }
    }
    return $biggest;
  }
  return;
}

=head2 get_display_xref

  DEPRECATED use $gene->display_xref

=cut

sub get_display_xref {
  my ($self, $gene) = @_;

  deprecate("display xref should retrieved from Gene object directly");

  if (!defined $gene) {
    throw("Must call with a Gene object");
  }

  my $sth = $self->prepare(
    qq(
      SELECT e.db_name,
             x.display_label,
             x.xref_id
      FROM   gene g, 
             xref x, 
             external_db e
      WHERE  g.gene_id = ?
        AND  g.display_xref_id = x.xref_id
        AND  x.external_db_id = e.external_db_id
  ));

  $sth->bind_param(1, $gene->dbID, SQL_INTEGER);
  $sth->execute();

  my ($db_name, $display_label, $xref_id) = $sth->fetchrow_array();
  if (!defined $xref_id) {
    return undef;
  }

  my $db_entry = Bio::EnsEMBL::DBEntry->new(-dbid       => $xref_id,
                                            -adaptor    => $self->db->get_DBEntryAdaptor(),
                                            -dbname     => $db_name,
                                            -display_id => $display_label);

  return $db_entry;
} ## end sub get_display_xref

=head2 get_description

  DEPRECATED, use gene->get_description

=cut

sub get_description {
  my ($self, $dbID) = @_;

  deprecate("Gene description should be loaded on gene retrieval. Use gene->get_description()");

  if (!defined $dbID) {
    throw("must call with dbID");
  }

  my $sth = $self->prepare(
    "SELECT description 
                            FROM   gene_description 
                            WHERE  gene_id = ?");

  $sth->bind_param(1, $dbID, SQL_INTEGER);
  $sth->execute();

  my @array = $sth->fetchrow_array();
  return $array[0];
}

=head2 fetch_by_Peptide_id

  DEPRECATED, use fetch_by_translation_stable_id()

=cut

sub fetch_by_Peptide_id {
  my ($self, $translation_stable_id) = @_;

  deprecate("Please use better named fetch_by_translation_stable_id \n" . caller(2));

  $self->fetch_by_translation_stable_id($translation_stable_id);
}

=head2 get_stable_entry_info

  DEPRECATED use $gene->stable_id instead

=cut

sub get_stable_entry_info {
  my ($self, $gene) = @_;

  deprecated("stable id info is loaded on default, no lazy loading necessary");

  if (!defined $gene || !ref $gene || !$gene->isa('Bio::EnsEMBL::Gene')) {
    throw("Needs a gene object, not a $gene");
  }

  my $created_date  = $self->db->dbc->from_date_to_seconds("created_date");
  my $modified_date = $self->db->dbc->from_date_to_seconds("modified_date");

  my $sth = $self->prepare("SELECT stable_id, " . $created_date . "," . $modified_date . ", version FROM gene WHERE gene_id = ?");

  $sth->bind_param(1, $gene->dbID, SQL_INTEGER);
  $sth->execute();

  my @array = $sth->fetchrow_array();
  $gene->{'stable_id'} = $array[0];
  $gene->{'created'}   = $array[1];
  $gene->{'modified'}  = $array[2];
  $gene->{'version'}   = $array[3];

  return 1;
} ## end sub get_stable_entry_info

=head2 fetch_all_by_DBEntry

  DEPRECATED - Use fetch_all_by_external_name instead

=cut

sub fetch_all_by_DBEntry {
  my $self = shift;

  deprecate('Use fetch_all_by_external_name instead.');

  return $self->fetch_all_by_external_name(@_);
}

*/

