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

#include "TranscriptAdaptor.h"
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
#include "MetaCoordContainer.h"
#include "BaseFeatureAdaptor.h"
#include "AttributeAdaptor.h"
#include "IntronSupportingEvidenceAdaptor.h"
#include "TranscriptSupportingFeatureAdaptor.h"
#include "TranslationAdaptor.h"

#include "ExonAdaptor.h"
#include "SliceAdaptor.h"

#include "StatementHandle.h"
#include "ResultRow.h"

#include "Error.h"
/*
=head1 DESCRIPTION

This adaptor provides a means to retrieve and store information related
to Transcripts.  Primarily this involves the retrieval or storage of
Bio::EnsEMBL::Transcript objects from a database.

See Bio::EnsEMBL::Transcript for details of the Transcript class.

=cut
*/

NameTableType TranscriptAdaptor_tableNames = {{"transcript","t"},
                                              {"xref","x"},
                                              {"external_db","exdb"},
                                              {NULL,NULL}};

NameTableType TranscriptAdaptor_leftJoins = {{"xref","x.xref_id = t.display_xref_id"},
                                             {"external_db","exdb.external_db_id = x.external_db_id"},
                                             {NULL,NULL}};




TranscriptAdaptor *TranscriptAdaptor_new(DBAdaptor *dba) {
  TranscriptAdaptor *ta;

  if ((ta = (TranscriptAdaptor *)calloc(1,sizeof(TranscriptAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for TranscriptAdaptor\n");
    return NULL;
  }
  BaseFeatureAdaptor_init((BaseFeatureAdaptor *)ta, dba, TRANSCRIPT_ADAPTOR);

  ta->getTables                  = TranscriptAdaptor_getTables;
  ta->getColumns                 = TranscriptAdaptor_getColumns;
  ta->store                      = (BaseAdaptor_StoreFunc)TranscriptAdaptor_store;
  ta->objectsFromStatementHandle = (BaseAdaptor_ObjectsFromStatementHandleFunc)TranscriptAdaptor_objectsFromStatementHandle;
  ta->leftJoin                   = TranscriptAdaptor_leftJoin;

  return ta;
}

/*
# _tables
#
#  Description: PROTECTED implementation of superclass abstract method.
#               Returns the names, aliases of the tables to use for queries.
#  Returntype : list of listrefs of strings
#  Exceptions : none
#  Caller     : internal
#  Status     : Stable
*/
NameTableType *TranscriptAdaptor_getTables() {
  return &TranscriptAdaptor_tableNames;
}


/*
#_columns
#
#  Description: PROTECTED implementation of superclass abstract method.
#               Returns a list of columns to use for queries.
#  Returntype : list of strings
#  Exceptions : none
#  Caller     : internal
#  Status     : Stable
*/
char *Transcript_cols[] = {
                           "t.transcript_id",
                           "t.seq_region_id",
                           "t.seq_region_start",
                           "t.seq_region_end",
                           "t.seq_region_strand",
                           "t.analysis_id",
                           "t.gene_id",
                           "t.is_current",
                           "t.stable_id",
                           "t.version",
                           "UNIX_TIMESTAMP(t.created_date)",
                           "UNIX_TIMESTAMP(t.modified_date)",
                           "t.description",
                           "t.biotype",
                           "exdb.db_name",
                           "exdb.status",
                           "exdb.db_display_name",
                           "x.xref_id",
                           "x.display_label",
                           "x.dbprimary_acc",
                           "x.version",
                           "x.description",
                           "x.info_type",
                           "x.info_text",
                           "exdb.db_release",
                           NULL };

char **TranscriptAdaptor_getColumns() {
  return Transcript_cols;
}

NameTableType *TranscriptAdaptor_leftJoin() {
  return &TranscriptAdaptor_leftJoins;
}

/*
=head2 fetch_by_stable_id

  Arg [1]    : String $stable_id 
               The stable id of the transcript to retrieve
  Example    : my $tr = $tr_adaptor->fetch_by_stable_id('ENST00000309301');
  Description: Retrieves a transcript via its stable id.
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/
Transcript *TranscriptAdaptor_fetchByStableId(TranscriptAdaptor *ta, char *stableId) {
  char constraint[1024];

  sprintf(constraint, "t.stable_id = '%s' AND t.is_current = 1", stableId);

  Vector *transcripts = TranscriptAdaptor_genericFetch(ta, constraint, NULL, NULL);
  Transcript *transcript = Vector_getElementAt(transcripts, 0);
  Vector_free(transcripts);

  return transcript;
}


Vector *TranscriptAdaptor_fetchAll(TranscriptAdaptor *ta) {
  char constraint[1024];

// NIY: Maybe just a constant string
  sprintf(constraint, "t.biotype != 'LRG_gene' and t.is_current = 1");

  return TranscriptAdaptor_genericFetch(ta, constraint, NULL, NULL);
}

/*
=head2 fetch_all_versions_by_stable_id 

  Arg [1]     : String $stable_id 
                The stable ID of the transcript to retrieve
  Example     : my $tr = $tr_adaptor->fetch_all_version_by_stable_id
                  ('ENST00000309301');
  Description : Similar to fetch_by_stable_id, but retrieves all versions of a
                transcript stored in the database.
  Returntype  : listref of Bio::EnsEMBL::Transcript objects
  Exceptions  : if we cant get the gene in given coord system
  Caller      : general
  Status      : At Risk

=cut
*/
Vector *TranscriptAdaptor_fetchAllVersionsByStableId(TranscriptAdaptor *ta, char *stableId) {
  char constraint[1024];

  sprintf(constraint, "t.stable_id = '%s'", stableId);

  return TranscriptAdaptor_genericFetch(ta, constraint, NULL, NULL);
}

/*
=head2 fetch_by_translation_stable_id

  Arg [1]    : String $transl_stable_id
               The stable identifier of the translation of the transcript to 
               retrieve
  Example    : my $tr = $tr_adaptor->fetch_by_translation_stable_id
                  ('ENSP00000311007');
  Description: Retrieves a Transcript object using the stable identifier of
               its translation.
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/
Transcript *TranscriptAdaptor_fetchByTranslationStableId(TranscriptAdaptor *ta, char *translationStableId) {
  char qStr[1024];

  sprintf(qStr, "SELECT  tr.transcript_id "
                "FROM    transcript tr, translation tl "
                "WHERE   tl.stable_id = '%s' "
                "AND     tr.transcript_id = tl.transcript_id "
                "AND     tr.is_current = 1", translationStableId);

  StatementHandle *sth = ta->prepare((BaseAdaptor *)ta,qStr,strlen(qStr));
  sth->execute(sth);

  if (sth->numRows(sth) == 0) {
    return NULL;
  }

  ResultRow *row = sth->fetchRow(sth);
  IDType transcriptId = row->getLongLongAt(row, 0);

  sth->finish(sth);

  Transcript *transcript = (Transcript *)TranscriptAdaptor_fetchByDbID(ta, transcriptId);
  return transcript;
}

/*
=head2 fetch_by_translation_id

  Arg [1]    : Int $id
               The internal identifier of the translation whose transcript
               is to be retrieved
  Example    : my $tr = $tr_adaptor->fetch_by_translation_id($transl->dbID);
  Description: Given the internal identifier of a translation this method 
               retrieves the transcript associated with that translation.
               If the transcript cannot be found undef is returned instead.
  Returntype : Bio::EnsEMBL::Transcript or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/
Transcript *TranscriptAdaptor_fetchByTranslationId(TranscriptAdaptor *ta, IDType translationId) {
  char qStr[1024];

  sprintf(qStr, "SELECT  transcript_id "
                "FROM    translation  "
                "WHERE   tl.translation_id = "IDFMTSTR,
                translationId);

  StatementHandle *sth = ta->prepare((BaseAdaptor *)ta,qStr,strlen(qStr));
  sth->execute(sth);

  if (sth->numRows(sth) == 0) {
    return NULL;
  }

  ResultRow *row = sth->fetchRow(sth);
  IDType transcriptId = row->getLongLongAt(row, 0);

  sth->finish(sth);

  Transcript *transcript = (Transcript *)TranscriptAdaptor_fetchByDbID(ta, transcriptId);
  return transcript;
}

/*
=head2 fetch_all_by_Gene

  Arg [1]    : Bio::EnsEMBL::Gene $gene
               The gene to fetch transcripts of
  Example    : my $gene = $gene_adaptor->fetch_by_stable_id('ENSG0000123');
               my @transcripts = { $tr_adaptor->fetch_all_by_Gene($gene) };
  Description: Retrieves Transcript objects for given gene. Puts Genes slice
               in each Transcript. 
  Returntype : Listref of Bio::EnsEMBL::Transcript objects
  Exceptions : none
  Caller     : Gene->get_all_Transcripts()
  Status     : Stable

=cut
*/
Vector *TranscriptAdaptor_fetchAllByGene(TranscriptAdaptor *ta, Gene *gene) {
  char constraint[1024];

  sprintf(constraint, "t.gene_id = "IDFMTSTR, Gene_getDbID(gene));

  // Use the fetch_all_by_Slice_constraint method because it handles the
  // difficult Haps/PARs and coordinate remapping.

  // Get a slice that entirely overlaps the gene.  This is because we
  // want all transcripts to be retrieved, not just ones overlapping
  // the slice the gene is on (the gene may only partially overlap the
  // slice).  For speed reasons, only use a different slice if necessary
  // though.

  Slice *gSlice = Gene_getSlice(gene);

  if (gSlice == NULL) {
    fprintf(stderr, "Gene must have attached slice to retrieve transcripts.\n");
    exit(1);
  }

  Slice *slice;

  if ( Gene_getStart(gene) < 1 || Gene_getEnd(gene) > Slice_getLength(gSlice) ) {
// No circular stuff
//    if ( $gslice->is_circular() ) {
//      $slice = $gslice;
//    } else {
    SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(ta->dba);
    slice = SliceAdaptor_fetchByFeature(sa, (SeqFeature *)gene, 0, 0);
//    }
  } else {
    slice = gSlice;
  }

  Vector *transcripts = TranscriptAdaptor_fetchAllBySliceConstraint(ta, slice, constraint, NULL);

  if ( slice != gSlice ) {
    Vector *out = Vector_new();
    int i;
    for (i=0; i<Vector_getNumElement(transcripts); i++) {
      Transcript *tr = Vector_getElementAt(transcripts, i);
      Vector_addElement(out, Transcript_transfer(tr, gSlice));
    }
// NIY Do I need to set a free func here - probably???
    Vector_free(transcripts);

    transcripts = out;
  }

  fprintf(stderr,"Note: Setting of canonical transcript not implemented yet in TranscriptAdaptor_fetchAllByGene\n");
// This is a nuts way of doing this!!!!
/*
  my $canonical_t = $gene->canonical_transcript();

  foreach my $t ( @{$transcripts} ) {
    if ( $t->equals($canonical_t) ) {
      $t->is_canonical(1);
      last;
    }
  }
*/
  return transcripts;
}


/*
=head2 fetch_all_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               The slice to fetch transcripts on
  Arg [2]    : (optional) Boolean $load_exons
               If true, exons will be loaded immediately rather than
               lazy loaded later
  Arg [3]    : (optional) String $logic_name
               The logic name of the type of features to obtain
  ARG [4]    : (optional) String $constraint
               An extra contraint.
  Example    : my @transcripts = @{ $tr_adaptor->fetch_all_by_Slice($slice) };
  Description: Overrides superclass method to optionally load exons
               immediately rather than lazy-loading them later. This
               is more efficient when there are a lot of transcripts whose
               exons are going to be used.
  Returntype : Listref of Bio::EnsEMBL::Transcript objects
  Exceptions : thrown if exon cannot be placed on transcript slice
  Caller     : Slice::get_all_Transcripts
  Status     : Stable

=cut
*/
typedef struct TranscriptRankPairStruct {
  Transcript *transcript;
  int rank;
} TranscriptRankPair;

TranscriptRankPair *TranscriptRankPair_new(Transcript *transcript, int rank) {
  TranscriptRankPair *trp;

  if ((trp = (TranscriptRankPair *)calloc(1,sizeof(TranscriptRankPair))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for TranscriptRankPair\n");
    return NULL;
  }
  trp->transcript = transcript;
  trp->rank = rank;

  return trp;
}

void TranscriptRankPair_free(TranscriptRankPair *trp) {
  free(trp);
}

Vector *TranscriptAdaptor_fetchAllBySlice(TranscriptAdaptor *ta, Slice *slice, int loadExons, char *logicName, char *inputConstraint) {
  char *constraint = NULL;
  if ((constraint = (char *)calloc(655500,sizeof(char))) == NULL) {
    fprintf(stderr,"Failed allocating constraint\n");
    return NULL;
  }

  strcpy(constraint, "t.is_current = 1");
  //fprintf(stderr, "Length of input constraint = %ld\n", strlen(inputConstraint));

  if (inputConstraint != NULL && inputConstraint[0] != '\0') {
    sprintf(constraint,"%s AND %s", constraint, inputConstraint);
  }
    
  Vector *transcripts = TranscriptAdaptor_fetchAllBySliceConstraint(ta, slice, constraint, logicName);

  // if there are 0 or 1 transcripts still do lazy-loading
// SMJS Tweaked so never does lazy loading
  if ( ! loadExons || Vector_getNumElement(transcripts) < 1 ) {
    return transcripts;
  }

  // preload all of the exons now, instead of lazy loading later
  // faster than 1 query per transcript

  // first check if the exons are already preloaded
  // @todo FIXME: Should test all exons.
  Transcript *firstTranscript = Vector_getElementAt(transcripts, 0);
  // Deliberate direct reach into transcript struct to avoid any lazy loading
  // Not sure if I should set exons to null when there aren't any loaded - maybe, but it makes everything harder
  if (Vector_getNumElement(firstTranscript->exons) != 0) {
    return transcripts;
  }

  // Get extent of region spanned by transcripts.
  long minStart =  2000000000;
  long maxEnd   = -2000000000;

  int i;
  for (i=0; i<Vector_getNumElement(transcripts); i++) {
    Transcript *t  = Vector_getElementAt(transcripts, i);
    if (Transcript_getSeqRegionStart((SeqFeature*)t) < minStart) {
      minStart = Transcript_getSeqRegionStart((SeqFeature*)t);
    }
    if (Transcript_getSeqRegionEnd((SeqFeature*)t) > maxEnd) {
      maxEnd = Transcript_getSeqRegionEnd((SeqFeature*)t);
    }
  }

  Slice *extSlice;

  if (minStart >= Slice_getStart(slice) && maxEnd <= Slice_getEnd(slice)) {
    extSlice = slice;
  } else {
    SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(ta->dba);
    extSlice = SliceAdaptor_fetchByRegion(sa, Slice_getCoordSystemName(slice), Slice_getSeqRegionName(slice),
                                          minStart, maxEnd, Slice_getStrand(slice), CoordSystem_getVersion(Slice_getCoordSystem(slice)), 0);
  }

  // associate exon identifiers with transcripts
  IDHash *trHash = IDHash_new(IDHASH_MEDIUM);
  for (i=0; i<Vector_getNumElement(transcripts); i++) {
    Transcript *t  = Vector_getElementAt(transcripts, i);
    if ( ! IDHash_contains(trHash, Transcript_getDbID(t))) {
      IDHash_add(trHash, Transcript_getDbID(t), t);
    }
  }

  IDType *uniqueIds = IDHash_getKeys(trHash);
  int nUniqueId = IDHash_getNumValues(trHash);

  int maxSize = 16384;

  char tmpStr[1024];
  char *qStr = NULL;
  int lenNum;
  IDHash *exTrHash = IDHash_new(IDHASH_LARGE);
  int endPoint;
//  bzero(qStr,655500);

  if ((qStr = (char *)calloc(655500,sizeof(char))) == NULL) {
    fprintf(stderr,"Failed allocating qStr\n");
    return transcripts;
  }


// Divide query if a lot of ids - Not done in perl
  for (i=0; i<nUniqueId; i+=maxSize) {
    //fprintf(stderr,"Transcript loop i = %d\n", i);
    endPoint = sprintf(qStr, "SELECT transcript_id, exon_id, rank FROM exon_transcript WHERE transcript_id IN (" );
    int j;
    for (j=0; j<maxSize && j+i<nUniqueId; j++) {
      if (j!=0) {
        qStr[endPoint++] = ',';
        qStr[endPoint++] = ' ';
      }
      lenNum = sprintf(tmpStr,IDFMTSTR,uniqueIds[i+j]);
      memcpy(&(qStr[endPoint]), tmpStr, lenNum);
      endPoint+=lenNum;
    }
    qStr[endPoint++] = ')';
    qStr[endPoint] = '\0';
  
    StatementHandle *sth = ta->prepare((BaseAdaptor *)ta,qStr,strlen(qStr));
    sth->execute(sth);
  
    ResultRow *row;
    while ((row = sth->fetchRow(sth))) {
      IDType trId = row->getLongLongAt(row,0);
      IDType exId = row->getLongLongAt(row,1);
      int    rank = row->getIntAt(row,2);
  
      if (! IDHash_contains(exTrHash, exId)) {
        Vector *vec = Vector_new();
        Vector_setFreeFunc(vec, TranscriptRankPair_free);
        IDHash_add(exTrHash, exId, vec);
      }
      Vector *exVec = IDHash_getValue(exTrHash, exId);
      TranscriptRankPair *trp = TranscriptRankPair_new(IDHash_getValue(trHash, trId), rank);
      Vector_addElement(exVec, trp);
    }

    sth->finish(sth);
  }

  free(uniqueIds);
  IDHash_free(trHash, NULL);

  //  sprintf( "e.exon_id IN (%s)",
  //    join( ',', sort { $a <=> $b } keys(%ex_tr_hash) ) ) );
  // Note this is a constraint rather than a complete query, but I'm using qStr to save stack space
  uniqueIds = IDHash_getKeys(exTrHash);
  nUniqueId = IDHash_getNumValues(exTrHash);

  qsort(uniqueIds, IDHash_getNumValues(exTrHash), sizeof(IDType), idTypeCompFunc); 

  Vector *exons = Vector_new();

//  bzero(qStr,655500);
  // Divide query if a lot of ids - Not done in perl
  for (i=0; i<nUniqueId; i+=maxSize) {
    //fprintf(stderr,"Exon loop i = %d\n", i);
    endPoint = sprintf(qStr, "e.exon_id IN (");
    int j;
    for (j=0; j<maxSize && j+i<nUniqueId; j++) {
      if (j!=0) {
        qStr[endPoint++] = ',';
        qStr[endPoint++] = ' ';
      }
      
      lenNum = sprintf(tmpStr,IDFMTSTR,uniqueIds[j+i]);
      memcpy(&(qStr[endPoint]), tmpStr, lenNum);
      endPoint+=lenNum;
    }
    qStr[endPoint++] = ')';
    qStr[endPoint] = '\0';

    //fprintf(stderr, "qStr = %s\n", qStr);
  
    // Interaction with slice feature fetch cache can be horrid - it frees the oldest cached features vector (and the features!) after cachce fills
   
    ExonAdaptor *ea = DBAdaptor_getExonAdaptor(ta->dba);
    Vector *tmpVec = ExonAdaptor_fetchAllBySliceConstraint(ea, extSlice, qStr, NULL);  
    
    //fprintf(stderr,"Adding %d elements from tmpVec to exons. Num in exons before = %d\n", Vector_getNumElement(tmpVec), Vector_getNumElement(exons));
    
    Vector_append(exons, tmpVec);
    //fprintf(stderr,"Num in exons after = %d\n", Vector_getNumElement(exons));

    // Interaction with slice feature fetch caching is nasty - don't try to free this vector
      //Vector_setFreeFunc(tmpVec, NULL);
      //Vector_free(tmpVec);
  }

    // move exons onto transcript slice, and add them to transcripts
  for (i=0; i<Vector_getNumElement(exons); i++) {
    Exon *ex = Vector_getElementAt(exons, i);
  
    // Perl didn't have this line - it was in GeneAdaptor version so I think I'm going to keep it
    if (!IDHash_contains(exTrHash, Exon_getDbID(ex))) {
      //fprintf(stderr,"Exon " IDFMTSTR " not found in exTrHash\n", Exon_getDbID(ex));
      continue;
    }
  
    Exon *newEx;
    if (slice != extSlice) {
      newEx = Exon_transfer(ex, slice);
      if (newEx == NULL) {
        fprintf(stderr, "Unexpected. Exon could not be transferred onto Transcript slice.\n");
        exit(1);
      }
    } else {
      newEx = ex;
    }
  
    Vector *exVec = IDHash_getValue(exTrHash, Exon_getDbID(newEx));
    int j;
    for (j=0; j<Vector_getNumElement(exVec); j++) {
      TranscriptRankPair *trp = Vector_getElementAt(exVec, j);
      //fprintf(stderr,"Adding exon "IDFMTSTR" at rank %d to transcript %p\n", Exon_getDbID(newEx), trp->rank, trp->transcript);
      Transcript_addExon(trp->transcript, newEx, trp->rank);
    }
  }

  Vector_free(exons);

  free(uniqueIds);

  TranslationAdaptor *tla = DBAdaptor_getTranslationAdaptor(ta->dba);

  // load all of the translations at once
  Vector *tmp = TranslationAdaptor_fetchAllByTranscriptList(tla, transcripts);

  // don't actually want this Vector so free (work was done in fetchAllByTranscriptList)
  Vector_free(tmp);
  //fprintf(stderr,"Translation fetching not yet implemented in transcript adaptor\n");

  // Free stuff
  IDHash_free(exTrHash, Vector_free);

  free(qStr);
  free(constraint);
  
  return transcripts;
}


/* Don't bother
=head2 fetch_all_by_external_name

  Arg [1]    : String $external_name
               An external identifier of the transcript to be obtained
  Arg [2]    : (optional) String $external_db_name
               The name of the external database from which the
               identifier originates.
  Arg [3]    : Boolean override. Force SQL regex matching for users
               who really do want to find all 'NM%'
  Example    : my @transcripts =
                  @{ $tr_adaptor->fetch_all_by_external_name( 'NP_065811.1') };
               my @more_transcripts = 
                  @{$tr_adaptor->fetch_all_by_external_name( 'NP_0658__._')};
  Description: Retrieves all transcripts which are associated with
               an external identifier such as a GO term, Swissprot
               identifer, etc.  Usually there will only be a single
               transcript returned in the list reference, but not
               always.  Transcripts are returned in their native
               coordinate system, i.e. the coordinate system in which
               they are stored in the database.  If they are required
               in another coordinate system the Transcript::transfer or
               Transcript::transform method can be used to convert them.
               If no transcripts with the external identifier are found,
               a reference to an empty list is returned.
               SQL wildcards % and _ are supported in the $external_name
               but their use is somewhat restricted for performance reasons.
               Users that really do want % and _ in the first three characters
               should use argument 3 to prevent optimisations
  Returntype : listref of Bio::EnsEMBL::Transcript
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_external_name {
  my ( $self, $external_name, $external_db_name, $override) = @_;

  my $entryAdaptor = $self->db->get_DBEntryAdaptor();

  my @ids =
    $entryAdaptor->list_transcript_ids_by_extids( $external_name,
                                                  $external_db_name, $override );

  return $self->fetch_all_by_dbID_list( \@ids );
}

=head2 fetch_all_by_GOTerm

  Arg [1]   : Bio::EnsEMBL::OntologyTerm
              The GO term for which transcripts should be fetched.

  Example:  @transcripts = @{
              $transcript_adaptor->fetch_all_by_GOTerm(
                $go_adaptor->fetch_by_accession('GO:0030326') ) };

  Description   : Retrieves a list of transcripts that are
                  associated with the given GO term, or with any of
                  its descendent GO terms.  The transcripts returned
                  are in their native coordinate system, i.e. in
                  the coordinate system in which they are stored
                  in the database.  If another coordinate system
                  is required then the Transcript::transfer or
                  Transcript::transform method can be used.

  Return type   : listref of Bio::EnsEMBL::Transcript
  Exceptions    : Throws of argument is not a GO term
  Caller        : general
  Status        : Stable

=cut

sub fetch_all_by_GOTerm {
  my ( $self, $term ) = @_;

  assert_ref( $term, 'Bio::EnsEMBL::OntologyTerm' );
  if ( $term->ontology() ne 'GO' ) {
    throw('Argument is not a GO term');
  }

  my $entryAdaptor = $self->db->get_DBEntryAdaptor();

  my %unique_dbIDs;
  foreach my $accession ( map { $_->accession() }
                          ( $term, @{ $term->descendants() } ) )
  {
    my @ids =
      $entryAdaptor->list_transcript_ids_by_extids( $accession, 'GO' );
    foreach my $dbID (@ids) { $unique_dbIDs{$dbID} = 1 }
  }

  my @result = @{
    $self->fetch_all_by_dbID_list(
                              [ sort { $a <=> $b } keys(%unique_dbIDs) ]
    ) };

  return \@result;
} ## end sub fetch_all_by_GOTerm

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
  my ( $self, $accession ) = @_;

  if ( $accession !~ /^GO:/ ) {
    throw('Argument is not a GO term accession');
  }

  my $goAdaptor =
    Bio::EnsEMBL::Registry->get_adaptor( 'Multi', 'Ontology',
                                         'OntologyTerm' );

  my $term = $goAdaptor->fetch_by_accession($accession);

  return $self->fetch_all_by_GOTerm($term);
}
*/

/*
=head2 fetch_by_display_label

  Arg [1]    : String $label - display label of transcript to fetch
  Example    : my $tr = $tr_adaptor->fetch_by_display_label("BRCA2");
  Description: Returns the transcript which has the given display label or
               undef if there is none. If there are more than 1, only the first
               is reported.
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/
Transcript *TranscriptAdaptor_fetchByDisplayLabel(TranscriptAdaptor *ta, char *label) {
  char constraint[1024];

  sprintf(constraint, "x.display_label = '%s' AND t.is_current = 1", label);

  Vector *transcripts = TranscriptAdaptor_genericFetch(ta, constraint, NULL, NULL);
  Transcript *transcript = Vector_getElementAt(transcripts, 0);
  Vector_free(transcripts);

  return transcript;
}


/*
=head2 fetch_all_by_exon_stable_id

  Arg [1]    : String $stable_id 
               The stable id of an exon in a transcript
  Example    : my $tr = $tr_adaptor->fetch_all_by_exon_stable_id
                  ('ENSE00000309301');
  Description: Retrieves a list of transcripts via an exon stable id.
  Returntype : Listref of Bio::EnsEMBL::Transcript objects
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/
Vector *TranscriptAdaptor_fetchByExonStableId(TranscriptAdaptor *ta, char *stableId) {
  char qStr[1024];

  sprintf(qStr,
      "SELECT t.transcript_id "
        "FROM transcript as t, "
             "exon_transcript as et, "
             "exon as e "
       "WHERE t.transcript_id = et.transcript_id "
         "AND et.exon_id = e.exon_id "
         "AND e.stable_id = '%s' "
         "AND e.is_current = 1", stableId);

  StatementHandle *sth = ta->prepare((BaseAdaptor *)ta,qStr,strlen(qStr));
  sth->execute(sth);

  if (sth->numRows(sth) == 0) {
    fprintf(stderr, "Failed fetching transcript using exon stable id %s - returning NULL\n",  stableId);
    return NULL;
  }

  Vector *transVec = Vector_new();

  ResultRow *row;
  while ((row = sth->fetchRow(sth))) {
    IDType dbID = row->getLongLongAt(row, 0);
    Transcript *transcript = (Transcript *)TranscriptAdaptor_fetchByDbID(ta, dbID);
    if (transcript != NULL) Vector_addElement(transVec, transcript);
  }
  sth->finish(sth);

  return transVec;
}

/*
=head2 fetch_all_by_biotype 

  Arg [1]    : String $biotype 
               listref of $biotypes
               The biotype of the gene to retrieve. You can also have a reference
               to a list of biotypes in the event of needing several.
  Example    : $transcript = $transcript_adaptor->fetch_all_by_biotype('pseudogene'); 
               $transcript = $transcript_adaptor->fetch_all_by_biotype(['protein_coding','ambiguous_orf']);
  Description: Retrieves an array reference of transcript objects from the 
               database via its biotype or biotypes.
               The transcript will be retrieved in its native coordinate system
               (i.e. in the coordinate system it is stored in the database). 
               It may be converted to a different coordinate system through a 
               call to transform() or transfer(). If the transcript is not found
               undef is returned instead.
  Returntype : listref of Bio::EnsEMBL::Transcript
  Exceptions : if we cant get the transcript in given coord system
  Caller     : general
  Status     : Stable

=cut
*/
Vector *TranscriptAdaptor_fetchAllByBiotype(TranscriptAdaptor *ta, Vector *biotypes) {
  char constraint[1024];

  TranscriptAdaptor_biotypeConstraint(ta, biotypes, constraint);

  return TranscriptAdaptor_genericFetch(ta, constraint, NULL, NULL);
}



void TranscriptAdaptor_biotypeConstraint(TranscriptAdaptor *ta, Vector *biotypes, char *constraint) {
  if (biotypes == NULL || Vector_getNumElement(biotypes) == 0) {
    fprintf(stderr,"list of biotypes expected\n");
    exit(1);
  }

  if (Vector_getNumElement(biotypes) > 1) {
    strcpy(constraint, "t.biotype IN (");

    int i;
    for (i=0;i<Vector_getNumElement(biotypes); i++) {
      char *biotype = Vector_getElementAt(biotypes, i);

      if (i>0) {
        strcat(constraint, ", ");
      }
      sprintf(constraint, "%s'%s'", constraint, biotype);
    }
    strcat(constraint, ") and t.is_current = 1");

  } else { // just one
    char *biotype = Vector_getElementAt(biotypes, 0);

    sprintf(constraint, "t.biotype = '%s' and t.is_current = 1", biotype);
  }

  return;
}

/* NIY
=head2 store

  Arg [1]    : Bio::EnsEMBL::Transcript $transcript
               The transcript to be written to the database
  Arg [2]    : Int $gene_dbID
               The identifier of the gene that this transcript is associated 
               with
  Arg [3]    : DEPRECATED (optional) Int $analysis_id
               The analysis_id to use when storing this gene. This is for 
               backward compatibility only and used to fall back to the gene
               analysis_id if no analysis object is attached to the transcript
               (which you should do for new code).
  Example    : $transID = $tr_adaptor->store($transcript, $gene->dbID);
  Description: Stores a transcript in the database and returns the new
               internal identifier for the stored transcript.
  Returntype : Int 
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/

static int doneAltTransWarn = 0;
IDType TranscriptAdaptor_store(TranscriptAdaptor *ta, Transcript *transcript, IDType geneDbID, IDType analysisId) {

  if (transcript == NULL) {
    fprintf(stderr, "feature is NULL in Transcript_store\n");
    exit(1);
  }

  Class_assertType(CLASS_TRANSCRIPT, transcript->objectType);

  DBAdaptor *db = ta->dba;
  AnalysisAdaptor *analysisAdaptor = DBAdaptor_getAnalysisAdaptor(db);

  if (Transcript_isStored(transcript, db)) {
    fprintf(stderr, "Transcript ["IDFMTSTR"] is already stored in this database.\n", Transcript_getDbID(transcript) );
    return Transcript_getDbID(transcript);
  }

  // Force lazy-loading of exons and ensure coords are correct.
/* NIY!!!!!!!!!!!!!!
  $transcript->recalculate_coordinates();
*/

/* Default set to 1 in Transcript_new so no need for this check
  my $is_current = ( defined( $transcript->is_current() )
                     ? $transcript->is_current()
                     : 1 );
*/
  int isCurrent = Transcript_getIsCurrent(transcript);

  // store analysis
  Analysis *analysis = Transcript_getAnalysis(transcript);
  IDType newAnalysisId;

  if (analysis != NULL) {
    if (!Analysis_isStored(analysis, db)) {
      AnalysisAdaptor_store(analysisAdaptor, analysis);
    }
    newAnalysisId = Analysis_getDbID(analysis);

/* Think above is equivalent to this
    if ( $analysis->is_stored($db) ) {
      $new_analysis_id = $analysis->dbID;
    } else {
      $new_analysis_id = $db->get_AnalysisAdaptor->store($analysis);
    }
*/
  } else if (analysisId) {
    // Fall back to analysis passed in (usually from gene) if analysis
    // wasn't set explicitly for the transcript. This is deprectated
    // though.
    fprintf(stderr,"You should explicitly attach an analysis object to the Transcript.\n"
                   "Will fall back to Gene analysis, but this behaviour is deprecated.\n" );
    newAnalysisId = analysisId;
  } else {
    fprintf(stderr,"Need an analysis_id to store the Transcript.\n");
    exit(1);
  }

  //
  // Store exons - this needs to be done before the possible transfer
  // of the transcript to another slice (in _prestore()).  Transfering
  // results in copies being made of the exons and we need to preserve
  // the object identity of the exons so that they are not stored twice
  // by different transcripts.
  // 
  ExonAdaptor *exonAdaptor = DBAdaptor_getExonAdaptor(db);
  int i;
  for (i=0; i<Transcript_getExonCount(transcript); i++) {
    Exon *exon = Transcript_getExonAt(transcript, i);
    ExonAdaptor_store(exonAdaptor, exon);
  }

/* Not doing transfers - don't think they should be necessary (hopefully)
  my $original_translation = $transcript->translation();
  my $original             = $transcript;
  my $seq_region_id;
  ( $transcript, $seq_region_id ) = $self->_pre_store($transcript);
*/
  IDType seqRegionId = BaseFeatureAdaptor_preStore((BaseFeatureAdaptor *)ta, (SeqFeature*)transcript);

  // First store the transcript without a display xref.  The display xref
  // needs to be set after xrefs are stored which needs to happen after
  // transcript is stored.

  //
  // Store transcript
  // 
  char fmtStr[1024];
  char qStr[1024];
  sprintf(fmtStr, "INSERT INTO transcript "
                "SET gene_id = "IDFMTSTR", "
                    "analysis_id = "IDFMTSTR", "
                    "seq_region_id = "IDFMTSTR", " 
                    "seq_region_start = %ld, "
                    "seq_region_end = %ld, "
                    "seq_region_strand = %d, " 
                    "biotype = %%s, "
                    "description = %%s, "
                    //"biotype = '%s', "
                    //"status = '%s', "
                    //"description = '%s', "
                    "is_current = %d, "
                    "canonical_translation_id = NULL",
         geneDbID,
         newAnalysisId,
         seqRegionId,
         Transcript_getSeqRegionStart((SeqFeature*)transcript),
          Transcript_getSeqRegionEnd((SeqFeature*)transcript),
          Transcript_getSeqRegionStrand((SeqFeature*)transcript),
         //Transcript_getBiotype(transcript),
         //Transcript_getStatus(transcript),
         //Transcript_getDescription(transcript),
         isCurrent); 

  char bioTypeQStr[1024];
  if (Transcript_getBiotype(transcript)) {
    sprintf(bioTypeQStr,"'%s'", Transcript_getBiotype(transcript));
  } else {
    sprintf(bioTypeQStr, "NULL");
  }
  char descQStr[1024];
  if (Transcript_getDescription(transcript)) {
    sprintf(descQStr,"'%s'", Transcript_getDescription(transcript));
  } else {
    sprintf(descQStr, "NULL");
  }

  sprintf(qStr, fmtStr, bioTypeQStr, descQStr);

  if (Transcript_getStableId(transcript)) {
/* Use FROM_UNIXTIME for now
    my $created  = $self->db->dbc->from_seconds_to_date($transcript->created_date());
    my $modified = $self->db->dbc->from_seconds_to_date($transcript->modified_date());
*/

    // Assume version will be positive, Transcript sets it to -1 when initialised
    int version = Transcript_getVersion(transcript) > 0 ? Transcript_getVersion(transcript) : 1; 
    sprintf(qStr,"%s, stable_id = '%s', version = %d, created_date = FROM_UNIXTIME(%ld), modified_date = FROM_UNIXTIME(%ld)",
            qStr, Transcript_getStableId(transcript), version, Transcript_getCreated(transcript), Transcript_getModified(transcript));
  }

  StatementHandle *tst = ta->prepare((BaseAdaptor *)ta,qStr,strlen(qStr));

  tst->execute(tst);

// Swapped id fetch and finish statements
  IDType transcDbID = tst->getInsertId(tst);

  tst->finish(tst);
  // I'll let the connection to finish before killing the Adaptor if the table_id is 0
  if (!transcDbID) {
    exit(1);
  }
  // 
  // Store translation
  //

  Vector *altTranslations = NULL;
/* NIY
  Vector *altTranslations = Transcript_getAllAlternativeTranslations(transcript);
*/

  Translation *translation = Transcript_getTranslation(transcript);

  if ( translation != NULL ) {
    // Make sure that the start and end exon are set correctly.
    Exon *startExon = Translation_getStartExon(translation);
    Exon *endExon   = Translation_getEndExon(translation);

    if ( startExon == NULL) {
      fprintf(stderr, "Translation does not define a start exon.\n");
      exit(1);
    }

    if ( endExon == NULL ) {
      fprintf(stderr, "Translation does not define an end exon.\n");
      exit(1);
    }

    // If the dbID is not set, this means the exon must have been a
    // different object in memory than the the exons of the transcript.
    // Try to find the matching exon in all of the exons we just stored.
// Use isStored instead
    if ( ! Exon_isStored(startExon, db) ) {
// NIY: Probably want to do comparisons rather than making string keys
      char key[2048];
      key[0] = '\0';
      Exon_getHashKey(startExon, key);

      Exon *newStartExon = NULL;
      for (i=0; i<Transcript_getExonCount(transcript) && !newStartExon; i++) {
        Exon *transExon = Transcript_getExonAt(transcript, i);
        char transExonKey[2048];
        Exon_getHashKey(transExon, transExonKey);
   
        if (!strcmp(transExonKey, key)) {
          newStartExon = transExon;
        }
      }
        
      if ( newStartExon ) {
        Translation_setStartExon(translation, newStartExon);
      } else {
        fprintf(stderr, "Translation's start_Exon does not appear "
                        "to be one of the exons in its associated Transcript.\n" );
        exit(1);
      }
    }

    if ( ! Exon_isStored(endExon, db) ) {
      char key[2048];
      key[0] = '\0';
      Exon_getHashKey(startExon, key);

      Exon *newEndExon = NULL;
      for (i=0; i<Transcript_getExonCount(transcript) && !newEndExon; i++) {
        Exon *transExon = Transcript_getExonAt(transcript, i);
        char transExonKey[2048];
        Exon_getHashKey(transExon, transExonKey);
   
        if (!strcmp(transExonKey, key)) {
          newEndExon = transExon;
        }
      }
        
      if ( newEndExon ) {
        Translation_setEndExon(translation, newEndExon);
      } else {
        fprintf(stderr, "Translation's end_Exon does not appear "
                        "to be one of the exons in its associated Transcript.\n");
        exit(1);
      }
    }

// Doesn't seem to be used    my $old_dbid = $translation->dbID();

    TranslationAdaptor *tlna = DBAdaptor_getTranslationAdaptor(db);
    TranslationAdaptor_store(tlna,  translation, transcDbID );

    // Need to update the canonical_translation_id for this transcript.

    sprintf(qStr, "UPDATE transcript SET canonical_translation_id = "IDFMTSTR" WHERE transcript_id = "IDFMTSTR,
            Translation_getDbID(translation), transcDbID);

    StatementHandle *sth = ta->prepare((BaseAdaptor *)ta,qStr,strlen(qStr));
    sth->execute(sth);
    sth->finish(sth);

    // Set values of the original translation, we may have copied it when
    // we transformed the transcript.
/* Shouldn't need because no transfer being done
    $original_translation->dbID( $translation->dbID() );
    $original_translation->adaptor( $translation->adaptor() );
*/
  }

  //
  // Store the alternative translations, if there are any.
  //

  if (!doneAltTransWarn) {
    fprintf(stderr, "Alt translation storing not implemented\n");
    doneAltTransWarn = 1;
  }
  if ( altTranslations != NULL &&
       Vector_getNumElement(altTranslations) > 0) {
    fprintf(stderr, "Alt translation storing not implemented\n");
/* NIY
    foreach my $alt_translation ( @{$alt_translations} ) {
      my $start_exon = $alt_translation->start_Exon();
      my $end_exon   = $alt_translation->end_Exon();

      if ( !defined($start_exon) ) {
        throw("Translation does not define a start exon.");
      } elsif ( !defined($end_exon) ) {
        throw("Translation does not defined an end exon.");
      }

      if ( !defined( $start_exon->dbID() ) ) {
        my $key = $start_exon->hashkey();
        ($start_exon) = grep { $_->hashkey() eq $key } @{$exons};

        if ( defined($start_exon) ) {
          $alt_translation->start_Exon($start_exon);
        } else {
          throw(   "Translation's start_Exon does not appear "
                 . "to be one of the exon in"
                 . "its associated Transcript" );
        }
      }
      if ( !defined( $end_exon->dbID() ) ) {
        my $key = $end_exon->hashkey();
        ($end_exon) = grep { $_->hashkey() eq $key } @$exons;

        if ( defined($end_exon) ) {
          $alt_translation->end_Exon($end_exon);
        } else {
          throw(   "Translation's end_Exon does not appear "
                 . "to be one of the exons in "
                 . "its associated Transcript." );
        }
      }

      $db->get_TranslationAdaptor()
        ->store( $alt_translation, $transc_dbID );
    }
*/
  }

  //
  // Store the xrefs/object xref mapping.
  // 
/* NIY
  DBEntryAdaptor *dbEntryAdaptor = DBAdaptor_getDBEntryAdaptor(db);

  Vector *dbEntries = Transcript_getAllDBEntries(transcript);
  for (i=0; i<Vector_getNumElement(dbEntries); i++) {
    DBEntry *dbe = Vector_getElementAt(dbEntries, i);
    DBEntryAdaptor_store(dbEntryAdaptor, dbe, transcDbID, "Transcript", 1);
  }


  //
  // Update transcript to point to display xref if it is set.
  //
  DBEntry *displayXref = Transcript_getDisplayXref(transcript);
  if (displayXref != NULL) {

    IDType dxrefId = 0;
    if (DBEntry_isStored(displayXref, db)) {
      dxrefId = DBEntry_getDbID(displayXref);
    } else {
      dxrefId = DBEntryAdaptor_exists(dbEntryAdaptor, displayXref);
    }

//    if (defined($dxref_id)) {
    if (dxrefId) {
      sprintf(qStr, "UPDATE transcript SET display_xref_id = "IDFMTSTR" WHERE transcript_id = "IDFMTSTR, dxrefId, transcDbID);

      StatementHandle *sth = ta->prepare((BaseAdaptor *)ta,qStr,strlen(qStr));

      sth->execute(sth);
      sth->finish(sth);

      DBEntry_setDbID(displayXref, dxrefId);
      DBEntry_setAdaptor(displayXref, (BaseAdaptor *)dbEntryAdaptor);
    } else {
      fprintf(stderr, "Display_xref %s:%s is not stored in database.\n"
                      "Not storing relationship to this transcript.\n",
              DBEntry_getDbName(displayXref), DBEntry_getDisplayId(displayXref));
      DBEntry_setDbID(displayXref, 0);
      DBEntry_setAdaptor(displayXref, NULL);
    }
  }
*/

  //
  // Link transcript to exons in exon_transcript table
  //
  sprintf(qStr, "INSERT INTO exon_transcript (exon_id,transcript_id,rank) " 
                "VALUES (%"IDFMTSTR",%"IDFMTSTR",%%d)" );

  StatementHandle *etst = ta->prepare((BaseAdaptor *)ta,qStr,strlen(qStr));

  int rank = 1;
  for (i=0; i<Transcript_getExonCount(transcript); i++) {
    Exon *exon = Transcript_getExonAt(transcript, i);

    etst->execute(etst, Exon_getDbID(exon), transcDbID, rank);
    rank++;
  }

  etst->finish(etst);

  // Now the supporting evidence
  TranscriptSupportingFeatureAdaptor *tsfAdaptor = DBAdaptor_getTranscriptSupportingFeatureAdaptor(db);
  TranscriptSupportingFeatureAdaptor_store(tsfAdaptor, transcDbID, Transcript_getAllSupportingFeatures(transcript));

  // store transcript attributes if there are any
  AttributeAdaptor *attrAdaptor = DBAdaptor_getAttributeAdaptor(db);
  Vector *attribs = Transcript_getAllAttributes(transcript, NULL);
  AttributeAdaptor_storeOnTranscriptId(attrAdaptor, transcDbID, attribs);
  Vector_free(attribs);

  // store the IntronSupportingEvidence features
  IntronSupportingEvidenceAdaptor *iseAdaptor = DBAdaptor_getIntronSupportingEvidenceAdaptor(db);
  Vector *intronSupportingEvidence = Transcript_getAllIntronSupportingEvidence(transcript);

  if (intronSupportingEvidence != NULL) {
    for (i=0; i<Vector_getNumElement(intronSupportingEvidence); i++) {
      IntronSupportingEvidence *ise = Vector_getElementAt(intronSupportingEvidence, i);
      IntronSupportingEvidenceAdaptor_store(iseAdaptor, ise);
      IntronSupportingEvidenceAdaptor_storeTranscriptLinkage(iseAdaptor, ise, transcript, transcDbID);
    }
  }

  // Update the original transcript object - not the transfered copy that
  // we might have created.
  Transcript_setDbID(transcript, transcDbID);
  Transcript_setAdaptor(transcript, (BaseAdaptor *)ta);

  return transcDbID;
}


/* Don't bother
=head2 get_Interpro_by_transid

  Arg [1]    : String $trans_stable_id
               The stable if of the transcript to obtain
  Example    : @i = $tr_adaptor->get_Interpro_by_transid($trans->stable_id()); 
  Description: Gets interpro accession numbers by transcript stable id.
               A hack really - we should have a much more structured 
               system than this.
  Returntype : listref of strings (Interpro_acc:description)
  Exceptions : none 
  Caller     : domainview? , GeneView
  Status     : Stable

=cut

sub get_Interpro_by_transid {
   my ($self,$trans_stable_id) = @_;

   my $sth = $self->prepare(qq(
      SELECT  STRAIGHT_JOIN i.interpro_ac, x.description
      FROM    transcript t,
              translation tl,
              protein_feature pf,
	      interpro i,
              xref x
      WHERE   t.stable_id = ?
      AND     tl.transcript_id = t.transcript_id
      AND     tl.translation_id = pf.translation_id
      AND     i.id = pf.hit_name
      AND     i.interpro_ac = x.dbprimary_acc
      AND     t.is_current = 1
  ));

  $sth->bind_param(1, $trans_stable_id, SQL_VARCHAR);
  $sth->execute();

  my @out;
  my %h;
  while( (my $arr = $sth->fetchrow_arrayref()) ) {
     if( $h{$arr->[0]} ) { next; }
     $h{$arr->[0]}=1;
     my $string = $arr->[0] .":".$arr->[1];
     push(@out,$string);
  }

  return \@out;
}
*/

/*
=head2 is_Transcript_canonical()

  Arg [1]     : Bio::EnsEMBL::Transcript $transcript
                The transcript to query with
  Example     : $tr_adaptor->is_Transcript_canonical($transcript);
  Description : Returns a boolean if the given transcript is considered
                canonical with respect to a gene
  Returntype  : Boolean
  Exceptions  : None
  Caller      : Bio::EnsEMBL::Transcript
  Status      : Beta
  

=cut
*/
int TranscriptAdaptor_isTranscriptCanonical(TranscriptAdaptor *ta, Transcript *transcript) {
  char qStr[1024];
  int flag = 0;

  sprintf(qStr, "select count(*) from gene where canonical_transcript_id = "IDFMTSTR, Transcript_getDbID(transcript));
  StatementHandle *sth = ta->prepare((BaseAdaptor *)ta,qStr,strlen(qStr));
  sth->execute(sth);

  if (sth->numRows(sth) > 0) {
    flag =1;
  }

  sth->finish(sth);
  return flag;
}


/* NIY
=head2 remove

  Arg [1]    : Bio::EnsEMBL::Transcript $transcript
               The transcript to remove from the database
  Example    : $tr_adaptor->remove($transcript);
  Description: Removes a transcript completely from the database, and all
               associated information.
               This method is usually called by the GeneAdaptor::remove method
               because this method will not preform the removal of genes
               which are associated with this transcript. Do not call this
               method directly unless you know there are no genes associated
               with the transcript!
  Returntype : none
  Exceptions : throw on incorrect arguments
               warning if transcript is not in this database
  Caller     : GeneAdaptor::remove
  Status     : Stable

=cut

sub remove {
  my $self = shift;
  my $transcript = shift;

  if(!ref($transcript) || !$transcript->isa('Bio::EnsEMBL::Transcript')) {
    throw("Bio::EnsEMBL::Transcript argument expected");
  }

  # sanity check: make sure nobody tries to slip past a prediction transcript
  # which inherits from transcript but actually uses different tables
  if($transcript->isa('Bio::EnsEMBL::PredictionTranscript')) {
    throw("TranscriptAdaptor can only remove Transcripts " .
          "not PredictionTranscripts");
  }

  if ( !$transcript->is_stored($self->db()) ) {
    warning("Cannot remove transcript ". $transcript->dbID .". Is not stored ".
            "in this database.");
    return;
  }

  # remove the supporting features of this transcript

  my $prot_adp = $self->db->get_ProteinAlignFeatureAdaptor;
  my $dna_adp = $self->db->get_DnaAlignFeatureAdaptor;

  my $sfsth = $self->prepare("SELECT feature_type, feature_id  " .
                             "FROM transcript_supporting_feature " .
                             "WHERE transcript_id = ?");

  $sfsth->bind_param(1, $transcript->dbID, SQL_INTEGER);
  $sfsth->execute();

  # statements to check for shared align_features
  my $sth1 = $self->prepare("SELECT count(*) FROM supporting_feature " .
			    "WHERE feature_type = ? AND feature_id = ?");
  my $sth2 = $self->prepare("SELECT count(*) " .
                            "FROM transcript_supporting_feature " .
			    "WHERE feature_type = ? AND feature_id = ?");

  SUPPORTING_FEATURE:
  while(my ($type, $feature_id) = $sfsth->fetchrow()){
    
    # only remove align_feature if this is the last reference to it
    $sth1->bind_param(1, $type, SQL_VARCHAR);
    $sth1->bind_param(2, $feature_id, SQL_INTEGER);
    $sth1->execute;
    $sth2->bind_param(1, $type, SQL_VARCHAR);
    $sth2->bind_param(2, $feature_id, SQL_INTEGER);
    $sth2->execute;
    my ($count1) = $sth1->fetchrow;
    my ($count2) = $sth2->fetchrow;
    if ($count1 + $count2 > 1) {
      #warn "transcript: shared feature, not removing $type|$feature_id\n";
      next SUPPORTING_FEATURE;
    }
    
    #warn "transcript: removing $type|$feature_id\n";
  
    if($type eq 'protein_align_feature'){
      my $f = $prot_adp->fetch_by_dbID($feature_id);
      $prot_adp->remove($f);
    }
    elsif($type eq 'dna_align_feature'){
      my $f = $dna_adp->fetch_by_dbID($feature_id);
      $dna_adp->remove($f);
    }
    else {
      warning("Unknown supporting feature type $type. Not removing feature.");
    }
  }
  $sfsth->finish();
  $sth1->finish();
  $sth2->finish();

  # delete the association to supporting features

  $sfsth = $self->prepare("DELETE FROM transcript_supporting_feature WHERE transcript_id = ?");
  $sfsth->bind_param(1, $transcript->dbID, SQL_INTEGER);
  $sfsth->execute();
  $sfsth->finish();
  
  # delete the associated IntronSupportingEvidence and if the ISE had no more
  # linked transcripts remove it
  my $ise_adaptor = $self->db->get_IntronSupportingEvidenceAdaptor();
  foreach my $ise (@{$transcript->get_all_IntronSupportingEvidence()}) {
    $ise_adaptor->remove_transcript_linkage($ise, $transcript);
    if(! $ise->has_linked_transcripts()) {
      $ise_adaptor->remove($ise);
    }
  }

  # remove all xref linkages to this transcript

  my $dbeAdaptor = $self->db->get_DBEntryAdaptor();
  foreach my $dbe (@{$transcript->get_all_DBEntries}) {
    $dbeAdaptor->remove_from_object($dbe, $transcript, 'Transcript');
  }

  # remove the attributes associated with this transcript
  my $attrib_adp = $self->db->get_AttributeAdaptor;  
  $attrib_adp->remove_from_Transcript($transcript);

  # remove the translation associated with this transcript

  my $translationAdaptor = $self->db->get_TranslationAdaptor();
  if( defined($transcript->translation()) ) {
    $translationAdaptor->remove( $transcript->translation );
  }

  # remove exon associations to this transcript

  my $exonAdaptor = $self->db->get_ExonAdaptor();
  foreach my $exon ( @{$transcript->get_all_Exons()} ) {
    # get the number of transcript references to this exon
    # only remove the exon if this is the last transcript to
    # reference it

    my $sth = $self->prepare( "SELECT count(*)
                               FROM   exon_transcript
                               WHERE  exon_id = ?" );
    $sth->bind_param(1, $exon->dbID, SQL_INTEGER);
    $sth->execute();
    my ($count) = $sth->fetchrow_array();
    $sth->finish();

    if($count == 1){
      $exonAdaptor->remove( $exon );
    }
  }

  my $sth = $self->prepare( "DELETE FROM exon_transcript
                             WHERE transcript_id = ?" );
  $sth->bind_param(1, $transcript->dbID, SQL_INTEGER);
  $sth->execute();
  $sth->finish();


  $sth = $self->prepare( "DELETE FROM transcript
                          WHERE transcript_id = ?" );
  $sth->bind_param(1, $transcript->dbID, SQL_INTEGER);
  $sth->execute();
  $sth->finish();

  $transcript->dbID(undef);
  $transcript->adaptor(undef);

  return;
}
*/


/* NIY
=head2 update

  Arg [1]    : Bio::EnsEMBL::Transcript $transcript
               The transcript to update
  Example    : $tr_adaptor->update($transcript);
  Description: Updates a transcript in the database.
  Returntype : None
  Exceptions : thrown if the $transcript is not a Bio::EnsEMBL::Transcript.
               warn if the method is called on a transcript that does not exist 
               in the database.
               Should warn if trying to update the number of attached exons, but
               this is a far more complex process and is not yet implemented.
  Caller     : general
  Status     : Stable

=cut

sub update {
  my ( $self, $transcript ) = @_;

  if (    !defined($transcript)
       || !ref($transcript)
       || !$transcript->isa('Bio::EnsEMBL::Transcript') )
  {
    throw("Must update a transcript object, not a $transcript");
  }

  my $update_transcript_sql = qq(
       UPDATE transcript
          SET analysis_id = ?,
              display_xref_id = ?,
              description = ?,
              biotype = ?,
              status = ?,
              is_current = ?,
              canonical_translation_id = ?
        WHERE transcript_id = ?
  );

  my $display_xref = $transcript->display_xref();
  my $display_xref_id;

  if ( defined($display_xref) && $display_xref->dbID() ) {
    $display_xref_id = $display_xref->dbID();
  } else {
    $display_xref_id = undef;
  }

  my $sth = $self->prepare($update_transcript_sql);

  $sth->bind_param( 1, $transcript->analysis()->dbID(), SQL_INTEGER );
  $sth->bind_param( 2, $display_xref_id, SQL_INTEGER );
  $sth->bind_param( 3, $transcript->description(), SQL_LONGVARCHAR );
  $sth->bind_param( 4, $transcript->biotype(),     SQL_VARCHAR );
  $sth->bind_param( 5, $transcript->status(),      SQL_VARCHAR );
  $sth->bind_param( 6, $transcript->is_current(),  SQL_TINYINT );
  $sth->bind_param( 7, (
                      defined( $transcript->translation() )
                      ? $transcript->translation()->dbID()
                      : undef ),
                    SQL_INTEGER );
  $sth->bind_param( 8, $transcript->dbID(), SQL_INTEGER );

  $sth->execute();
} ## end sub update
*/


/*
=head2 list_dbIDs

  Example    : @transcript_ids = @{ $t_adaptor->list_dbIDs };
  Description: Gets a list of internal ids for all transcripts in the db.
  Arg[1]     : <optional> int. not 0 for the ids to be sorted by the seq_region.  Returntype : Listref of Ints
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/
Vector *TranscriptAdaptor_listDbIDs(TranscriptAdaptor *ta, int ordered) {

// NIY: Shouldn't really do this direct BaseAdaptor call but ...
  return BaseAdaptor_listDbIDs((BaseAdaptor *)ta, "transcript", NULL, ordered);
}

/*
=head2 list_stable_ids

  Example    : @stable_trans_ids = @{ $transcript_adaptor->list_stable_ids };
  Description: Gets a list of stable ids for all transcripts in the current
               database.
  Returntype : Listref of Strings
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/
Vector *TranscriptAdaptor_listStableIDs(TranscriptAdaptor *ta) {

// NIY: Shouldn't really do this direct BaseAdaptor call but ...
  return BaseAdaptor_listDbIDs((BaseAdaptor *)ta, "transcript", "stable_id", 0);
}

/*
#_objs_from_sth

#  Arg [1]    : StatementHandle $sth
#  Arg [2]    : Bio::EnsEMBL::AssemblyMapper $mapper
#  Arg [3]    : Bio::EnsEMBL::Slice $dest_slice
#  Description: PROTECTED implementation of abstract superclass method.
#               Responsible for the creation of Transcripts.
#  Returntype : Listref of Bio::EnsEMBL::Transcripts in target coord system
#  Exceptions : none
#  Caller     : internal
#  Status     : Stable
*/
Vector *TranscriptAdaptor_objectsFromStatementHandle(TranscriptAdaptor *ta,
                                               StatementHandle *sth,
                                               AssemblyMapper *assMapper,
                                               Slice *destSlice) {
  //
  // This code is ugly because an attempt has been made to remove as many
  // function calls as possible for speed purposes.  Thus many caches and
  // a fair bit of gymnastics is used.
  //
  SliceAdaptor *sa     = DBAdaptor_getSliceAdaptor(ta->dba);
  AnalysisAdaptor *aa  = DBAdaptor_getAnalysisAdaptor(ta->dba);
  DBEntryAdaptor *dbea = DBAdaptor_getDBEntryAdaptor(ta->dba);

  Vector *transcripts = Vector_new();
  IDHash *sliceHash = IDHash_new(IDHASH_SMALL);

/* Don't bother with these three - analysis is cached in its adaptor, and I can't believe speed
   is going to be limited by name and cs access functions on a slice!
  my %analysis_hash;
  my %sr_name_hash;
  my %sr_cs_hash;
*/



/* Basically unused (is used but in an if else condition where both if and else do same thing)
   Not usre in GeneAdaptor version so think it can go
  my $asm_cs;
  my $cmp_cs;
  my $asm_cs_vers;
  my $asm_cs_name;
  my $cmp_cs_vers;
  my $cmp_cs_name;
  if($mapper) {
    $asm_cs = $mapper->assembled_CoordSystem();
    $cmp_cs = $mapper->component_CoordSystem();
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
  CoordSystem *destSliceCs;
  char * destSliceSrName;
  IDType destSliceSrId = 0;
  AssemblyMapperAdaptor *asma;

  if (destSlice) {
    destSliceStart  = Slice_getStart(destSlice);
    destSliceEnd    = Slice_getEnd(destSlice);
    destSliceStrand = Slice_getStrand(destSlice);
    destSliceLength = Slice_getLength(destSlice);
    destSliceCs     = Slice_getCoordSystem(destSlice);
    destSliceSrName = Slice_getSeqRegionName(destSlice);
    destSliceSrId   = Slice_getSeqRegionId(destSlice);
    asma            = DBAdaptor_getAssemblyMapperAdaptor(ta->dba);
  }


// Note FEATURE label is here
//FEATURE: while ($sth->fetch())
  ResultRow *row;
  while ((row = sth->fetchRow(sth))) {
    IDType transcriptId =    row->getLongLongAt(row,0);
    IDType seqRegionId =     row->getLongLongAt(row,1);
    long seqRegionStart =    row->getLongAt(row,2);
    long seqRegionEnd =      row->getLongAt(row,3);
    int seqRegionStrand =    row->getIntAt(row,4);
    IDType analysisId =      row->getLongLongAt(row,5);
    int isCurrent =          row->getIntAt(row,7);
    char *stableId =         row->getStringAt(row,8);
    int version =            row->getIntAt(row,9);
    int createdDate =        row->getIntAt(row,10);
    int modifiedDate =       row->getIntAt(row,11);
    char *description =      row->getStringAt(row,12);
    char *biotype =          row->getStringAt(row,13);
    char *externalDb =       row->getStringAt(row,14);
    char *externalStatus =   row->getStringAt(row,15);
    char *externalDbName =   row->getStringAt(row,16);
// Note changed from xrefId to displayXrefId to match GeneAdaptor version of this code
    IDType displayXrefId =   row->getLongLongAt(row,17);
    char *xrefDisplayLabel = row->getStringAt(row,18);
    char *xrefPrimaryAcc =   row->getStringAt(row,19);
    char *xrefVersion =      row->getStringAt(row,20);
// Note changed from xrefDescription to xrefDesc to match GeneAdaptor version of this code
    char *xrefDesc =         row->getStringAt(row,21);
    char *xrefInfoType =     row->getStringAt(row,22);
    char *xrefInfoText =     row->getStringAt(row,23);
// Added externalRelease to be consistent with Gene version
    char *externalRelease =  row->getStringAt(row,24);

    Analysis *analysis = AnalysisAdaptor_fetchByDbID(aa, analysisId);

    // Not doing internal seq id stuff for now
//    #need to get the internal_seq_region, if present
//    $seq_region_id = $self->get_seq_region_id_internal($seq_region_id);
//    my $slice = $slice_hash{"ID:".$seq_region_id};
//    my $dest_mapper = $mapper;
//
//    if(!$slice) {
//      $slice = $sa->fetch_by_seq_region_id($seq_region_id);
//      $slice_hash{"ID:".$seq_region_id} = $slice;
//      $sr_name_hash{$seq_region_id} = $slice->seq_region_name();
//      $sr_cs_hash{$seq_region_id} = $slice->coord_system();
//    }

    if (! IDHash_contains(sliceHash, seqRegionId)) {
      IDHash_add(sliceHash, seqRegionId, SliceAdaptor_fetchBySeqRegionId(sa, seqRegionId, POS_UNDEF, POS_UNDEF, STRAND_UNDEF));
    }
    Slice *slice = IDHash_getValue(sliceHash, seqRegionId);

    // obtain a mapper if none was defined, but a dest_seq_region was
    if (assMapper == NULL &&
        destSlice != NULL &&
        CoordSystem_compare(destSliceCs, Slice_getCoordSystem(slice))) {
      assMapper = AssemblyMapperAdaptor_fetchByCoordSystems(asma, destSliceCs, Slice_getCoordSystem(slice));
/*
      $asm_cs = $dest_mapper->assembled_CoordSystem();
      $cmp_cs = $dest_mapper->component_CoordSystem();
      $asm_cs_name = $asm_cs->name();
      $asm_cs_vers = $asm_cs->version();
      $cmp_cs_name = $cmp_cs->name();
      $cmp_cs_vers = $cmp_cs->version();
*/
    }

    Slice *transcriptSlice = slice;

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
      } else {
        mrs = AssemblyMapper_fastMap(assMapper, srName, seqRegionStart, seqRegionEnd, seqRegionStrand, srCs, NULL);
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


// was...
      //get a slice in the coord system we just mapped to
/* Identical code in if and else, so why test???
      if($asm_cs == $sr_cs || ($cmp_cs != $sr_cs && $asm_cs->equals($sr_cs))) {
        $slice = $slice_hash{"ID:".$seq_region_id} ||=
          $sa->fetch_by_seq_region_id($seq_region_id);
      } else {
        $slice = $slice_hash{"ID:".$seq_region_id} ||=
          $sa->fetch_by_seq_region_id($seq_region_id);
      }
*/
// Instead...
      if (! IDHash_contains(sliceHash, seqRegionId)) {
        IDHash_add(sliceHash, seqRegionId, SliceAdaptor_fetchBySeqRegionId(sa, seqRegionId, POS_UNDEF, POS_UNDEF, STRAND_UNDEF));
      }
      transcriptSlice = IDHash_getValue(sliceHash, seqRegionId);
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
/*
        if ( $dest_slice->is_circular ) 
          if ( $seq_region_start > $seq_region_end ) {
            # Looking at a feature overlapping the chromsome origin.
            if ( $seq_region_end > $dest_slice_start ) {
              # Looking at the region in the beginning of the chromosome
              $seq_region_start -= $dest_slice->seq_region_length();
            }
            if ( $seq_region_end < 0 ) {
              $seq_region_end += $dest_slice->seq_region_length();
            }
          } else {
            if (    $dest_slice_start > $dest_slice_end
                 && $seq_region_end < 0 )
            {
              # Looking at the region overlapping the chromosome
              # origin and a feature which is at the beginning of the
              # chromosome.
              $seq_region_start += $dest_slice->seq_region_length();
              $seq_region_end   += $dest_slice->seq_region_length();
            }
          }
*/
        }
      } else {
        // Negative strand.
// No circular stuff
        if (0) {
/*
        if (    $dest_slice->is_circular()
             && $seq_region_start > $seq_region_end )
          if ( $seq_region_end > $dest_slice_start ) {
            # Looking at the region in the beginning of the chromosome.
            $seq_region_start = $dest_slice_end - $seq_region_end + 1;
            $seq_region_end =
              $seq_region_end -
              $dest_slice->seq_region_length() -
              $dest_slice_start + 1;
          } else {
            my $tmp_seq_region_start = $seq_region_start;
            $seq_region_start =
              $dest_slice_end -
              $seq_region_end -
              $dest_slice->seq_region_length() + 1;
            $seq_region_end =
              $dest_slice_end - $tmp_seq_region_start + 1;
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
      transcriptSlice = destSlice;
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


    // Finally, create the new Transcript.
    Transcript *transcript = Transcript_new();
  
    Transcript_setAnalysis      (transcript, analysis);
    Transcript_setBiotype       (transcript, biotype);
    Transcript_setStart         (transcript, seqRegionStart);
    Transcript_setEnd           (transcript, seqRegionEnd);
    Transcript_setStrand        (transcript, seqRegionStrand);
    Transcript_setAdaptor       (transcript, (BaseAdaptor *)ta);
    Transcript_setSlice         (transcript, transcriptSlice);
    Transcript_setDbID          (transcript, transcriptId);
    Transcript_setStableId      (transcript, stableId);
    Transcript_setVersion       (transcript, version);
// Had $created_date || undef, for this (and equivalent for modified_date - not sure about that???
    Transcript_setCreated       (transcript, createdDate);
    Transcript_setModified      (transcript, modifiedDate);
    Transcript_setDescription   (transcript, description);
    Transcript_setExternalName  (transcript, xrefDisplayLabel); // Perl set this but shouldn't it be like gene with it coming from xref?
    Transcript_setExternalDb    (transcript, externalDb);
    Transcript_setExternalStatus(transcript, externalStatus);
    Transcript_setDisplayXref   (transcript, displayXref);
    Transcript_setIsCurrent     (transcript, isCurrent);
    Transcript_setEditsEnabled  (transcript, 1);

    Vector_addElement(transcripts, transcript);
    // Not sure what this one is????      'external_display_name' => $external_db_name,
    //  Doesn't seem to match anything in transcript object - there's external_db_name but that's not used in constructor
  }
  IDHash_free(sliceHash, NULL);

  return transcripts;
}


/* NIY
=head2 fetch_all_by_exon_supporting_evidence

  Arg [1]    : String $hit_name
               Name of supporting feature
  Arg [2]    : String $feature_type 
               one of "dna_align_feature" or "protein_align_feature"
  Arg [3]    : (optional) Bio::Ensembl::Analysis
  Example    : $tr = $tr_adaptor->fetch_all_by_exon_supporting_evidence
                  ('XYZ', 'dna_align_feature');
  Description: Gets all the transcripts with exons which have a specified hit
               on a particular type of feature. Optionally filter by analysis.
  Returntype : Listref of Bio::EnsEMBL::Transcript objects
  Exceptions : If feature_type is not of correct type.
  Caller     : general
  Status     : Stable

=cut
sub fetch_all_by_exon_supporting_evidence {
  my ($self, $hit_name, $feature_type, $analysis) = @_;

  if($feature_type !~ /(dna)|(protein)_align_feature/) {
    throw("feature type must be dna_align_feature or protein_align_feature");
  }

  my $anal_from = "";
  $anal_from = ", analysis a " if ($analysis);
  my $anal_where = "";
  $anal_where = "AND a.analysis_id = f.analysis_id AND a.analysis_id=? "
    if ($analysis);

  my $sql = qq(
      SELECT DISTINCT(t.transcript_id)
        FROM transcript t,
             exon_transcript et,
             supporting_feature sf,
             $feature_type f
             $anal_from
       WHERE t.transcript_id = et.transcript_id
         AND t.is_current = 1
         AND et.exon_id = sf.exon_id
         AND sf.feature_id = f.${feature_type}_id
         AND sf.feature_type = ?
         AND f.hit_name=?
         $anal_where
  );

  my $sth = $self->prepare($sql);

  $sth->bind_param(1, $feature_type, SQL_VARCHAR);
  $sth->bind_param(2, $hit_name, SQL_VARCHAR);
  $sth->bind_param(3, $analysis->dbID(), SQL_INTEGER) if ($analysis);

  $sth->execute();

  my @transcripts;

  while( my $id = $sth->fetchrow_array ) {
    my $transcript = $self->fetch_by_dbID( $id  );
    push(@transcripts, $transcript) if $transcript;
  }

  return \@transcripts;
}
*/


/* NIY
=head2 fetch_all_by_transcript_supporting_evidence

  Arg [1]    : String $hit_name
               Name of supporting feature
  Arg [2]    : String $feature_type 
               one of "dna_align_feature" or "protein_align_feature"
  Arg [3]    : (optional) Bio::Ensembl::Analysis
  Example    : $transcripts = $transcript_adaptor->fetch_all_by_transcript_supporting_evidence('XYZ', 'dna_align_feature');
  Description: Gets all the transcripts with evidence from a specified hit_name on a particular type of feature, stored in the  
               transcript_supporting_feature table. Optionally filter by analysis.  For hits stored in the supporting_feature 
               table (linked to exons) use fetch_all_by_exon_supporting_evidence instead.
  Returntype : Listref of Bio::EnsEMBL::Transcript objects
  Exceptions : If feature_type is not of correct type.
  Caller     : general
  Status     : Stable

=cut
sub fetch_all_by_transcript_supporting_evidence {
  
  my ($self, $hit_name, $feature_type, $analysis) = @_;

  if($feature_type !~ /(dna)|(protein)_align_feature/) {
    throw("feature type must be dna_align_feature or protein_align_feature");
  }

  my $anal_from = "";
  $anal_from = ", analysis a " if ($analysis);
  my $anal_where = "";
  $anal_where = "AND a.analysis_id = f.analysis_id AND a.analysis_id=? "
    if ($analysis);

  my $sql = qq(
      SELECT DISTINCT(t.transcript_id)
        FROM transcript t,
             transcript_supporting_feature sf,
             $feature_type f
             $anal_from
       WHERE t.transcript_id = sf.transcript_id
         AND t.is_current = 1
         AND sf.feature_id = f.${feature_type}_id
         AND sf.feature_type = ?
         AND f.hit_name=?
         $anal_where
  );

  my $sth = $self->prepare($sql);

  $sth->bind_param(1, $feature_type, SQL_VARCHAR);
  $sth->bind_param(2, $hit_name, SQL_VARCHAR);
  $sth->bind_param(3, $analysis->dbID(), SQL_INTEGER) if ($analysis);

  $sth->execute();

  my @transcripts;

  while( my $id = $sth->fetchrow_array ) {
    my $transcript = $self->fetch_by_dbID( $id  );
    push(@transcripts, $transcript) if $transcript;
  }

  return \@transcripts;
}
*/

/*

##########################
#                        #
#  DEPRECATED METHODS    #
#                        #
##########################


=head2 get_display_xref

  Description: DEPRECATED. Use $transcript->display_xref() instead.

=cut

sub get_display_xref {
  my ($self, $transcript) = @_;
	
  deprecate("display_xref should be retreived from Transcript object directly.");
  
  if ( !defined $transcript ) {
    throw("Must call with a Transcript object");
  }

  my $sth = $self->prepare(qq(
      SELECT e.db_name,
             x.display_label,
             e.db_external_name,
             x.xref_id
      FROM   transcript t, 
             xref x, 
             external_db e
      WHERE  t.transcript_id = ?
        AND  t.display_xref_id = x.xref_id
        AND  x.external_db_id = e.external_db_id
  ));
  
  $sth->bind_param(1, $transcript->dbID, SQL_INTEGER);
  $sth->execute();

  my ($db_name, $display_label, $xref_id, $display_db_name ) =
    $sth->fetchrow_array();
  
  if ( !defined $xref_id ) {
    return undef;
  }

  my $db_entry = Bio::EnsEMBL::DBEntry->new(
     -dbid => $xref_id,
     -adaptor => $self->db->get_DBEntryAdaptor(),
     -dbname => $db_name,
     -display_id => $display_label
     -db_display_name => $display_db_name
  );

  return $db_entry;
}


=head2 get_stable_entry_info

  Description: DEPRECATED. Use $transcript->stable_id() instead.

=cut

sub get_stable_entry_info {
  my ($self, $transcript) = @_;

  deprecate("Stable ids should be loaded directly now");

  unless ( defined $transcript && ref $transcript && 
	  $transcript->isa('Bio::EnsEMBL::Transcript') ) {
    throw("Needs a Transcript object, not a $transcript");
  }

  my $sth = $self->prepare(qq(
      SELECT stable_id, version 
      FROM   transcript
      WHERE  transcript_id = ?
  ));
                            
  $sth->bind_param(1, $transcript->dbID, SQL_INTEGER);
  $sth->execute();

  my @array = $sth->fetchrow_array();
  $transcript->{'_stable_id'} = $array[0];
  $transcript->{'_version'}   = $array[1];

  return 1;
}


=head2 fetch_all_by_DBEntry

  Description: DEPRECATED. Use fetch_all_by_external_name() instead.

=cut

sub fetch_all_by_DBEntry {
  my $self = shift;
  deprecate('Use fetch_all_by_external_name instead.');
  return $self->fetch_all_by_external_name(@_);
}
*/


