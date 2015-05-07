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

This adaptor provides a means to retrieve and store
Bio::EnsEMBL::Translation objects from/in a database.

Translation objects only truly make sense in the context of their
transcripts so the recommended means to retrieve Translations is
by retrieving the Transcript object first, and then fetching the
Translation.
*/
#include "TranslationAdaptor.h"
#include "TranscriptAdaptor.h"
#include "AttributeAdaptor.h"
#include "DBEntryAdaptor.h"
#include "BaseAdaptor.h"
#include "MysqlUtil.h"
#include "Exon.h"

#include "StatementHandle.h"
#include "ResultRow.h"


TranslationAdaptor *TranslationAdaptor_new(DBAdaptor *dba) {
  TranslationAdaptor *ta;

  if ((ta = (TranslationAdaptor *)calloc(1,sizeof(TranslationAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for TranslationAdaptor\n");
    exit(1);
  }
  BaseAdaptor_init((BaseAdaptor *)ta, dba, TRANSLATION_ADAPTOR);

  return ta;
}



/*
=head2 fetch_all_alternative_by_Transcript

  Arg [1]    : Bio::EnsEMBL::Transcript $transcript
  Example    :

    @tl = @{
      $translation_adaptor->fetch_all_alternative_by_Transcript(
                                                            $transcript)
      };

  Description: Retrieves all alternative translations associated with a
               particular transcript.  If no alternative translation is
               found, a reference to an empty list is returned.

  Returntype : listref of Bio::EnsEMBL::Translation
  Exceptions : throw on incorrect argument
  Caller     : Transcript
  Status     : Stable

=cut
*/
Vector *TranslationAdaptor_fetchAllAlternativeByTranscript(TranslationAdaptor *tlna, Transcript *transcript) {
  if (transcript == NULL) {
    fprintf(stderr, "NULL transcript in TranslationAdaptor_fetchByTranscript\n");
    exit(1);
  }

/*
  my $tl_created_date =
    $self->db()->dbc()->from_date_to_seconds('tl.created_date');
  my $tl_modified_date =
    $self->db()->dbc()->from_date_to_seconds('tl.modified_date');
*/
// Instead of the above...
  char *tlCreatedDate = "UNIX_TIMESTAMP(tl.createdDate)";
  char *tlModifiedDate = "UNIX_TIMESTAMP(tl.modifiedDate)";


  char qStr[1024];
  sprintf(qStr,"SELECT tl.translation_id, tl.start_exon_id, "
                      "tl.end_exon_id, tl.seq_start, tl.seq_end, "
                      "tl.stable_id, tl.version, %s, %s "
                 "FROM translation tl "
                 "JOIN transcript t "
                   "ON (t.transcript_id = tl.transcript_id) "
                "WHERE tl.transcript_id = "IDFMTSTR
                  "AND tl.translation_id != t.canonical_translation_id",
           tlCreatedDate, tlModifiedDate, Transcript_getDbID(transcript));

  StatementHandle *sth = tlna->prepare((BaseAdaptor *)tlna,qStr,strlen(qStr));

  sth->execute(sth);

  // Get all alternative translations.
  Vector *translations = Vector_new();

  ResultRow *row;
  while ((row = sth->fetchRow(sth))) {
    Translation *translation = TranslationAdaptor_translationFromResultRow(tlna, row, transcript);

// Huh, ???
//    $translation->transcript($transcript);

    Vector_addElement(translations, translation);

  }

  return translations;
}

/*
=head2 fetch_by_Transcript

  Arg [1]    : Bio::EnsEMBL::Transcript $transcript
  Example    : $tl = $translation_adaptor->fetch_by_Transcript($transcript);
  Description: Retrieves a Translation via its associated transcript.
               If the Translation is not found, undef is returned.
  Returntype : Bio::EnsEMBL::Translation
  Exceptions : throw on incorrect argument
  Caller     : Transcript
  Status     : Stable

=cut
*/
Translation *TranslationAdaptor_fetchByTranscript(TranslationAdaptor *tlna, Transcript *transcript) {

  if (transcript == NULL) {
    fprintf(stderr, "NULL transcript in TranslationAdaptor_fetchByTranscript\n");
    exit(1);
  }

/*
  my $tl_created_date =
    $self->db()->dbc()->from_date_to_seconds('tl.created_date');
  my $tl_modified_date =
    $self->db()->dbc()->from_date_to_seconds('tl.modified_date');
*/
// Instead of the above...
  char *tlCreatedDate = "UNIX_TIMESTAMP(tl.createdDate)";
  char *tlModifiedDate = "UNIX_TIMESTAMP(tl.modifiedDate)";

  char qStr[1024];
  sprintf(qStr, "SELECT tl.translation_id, tl.start_exon_id, "
                       "tl.end_exon_id, tl.seq_start, tl.seq_end, "
                       "tl.stable_id, tl.version, %s, %s "
                  "FROM translation tl "
                  "JOIN transcript tr "
                    "ON (tl.translation_id = tr.canonical_translation_id) "
                 "WHERE tr.transcript_id = "IDFMTSTR,
           tlCreatedDate, tlModifiedDate, Transcript_getDbID(transcript));

  StatementHandle *sth = tlna->prepare((BaseAdaptor *)tlna,qStr,strlen(qStr));

  sth->execute(sth);

  if (sth->numRows(sth) > 1) {
    fprintf(stderr, "Error: Expected one result when querying for translation with transcript\n");
    exit(1);
  } else if (sth->numRows(sth) == 0) {
    return NULL;
  }

  ResultRow *row = sth->fetchRow(sth);

  Translation *translation = TranslationAdaptor_translationFromResultRow(tlna, row, transcript); 

// Huh, don't think other fetch method did this??
//  $translation->transcript($transcript);

  return translation;
}

/* NIY
=head2 fetch_all_by_external_name

  Arg [1]    : string $external_name
               The external identifier for the translation(s) to be
               obtained.
  Arg [2]    : (optional) string $external_db_name
               The name of the external database from which the
               identifier originates.
  Arg [3]    : Boolean override. Force SQL regex matching for users
               who really do want to find all 'NM%'
  Example    : my @translations =
                  @{ $trl_adaptor->fetch_all_by_external_name('BRCA2') };
               my @many_translations = 
                  @{ $trl_adaptor->fetch_all_by_external_name('BRCA%') };
  Description: Retrieves a list of translations fetched via an
               external identifier.  Note that this may not be a
               particularly useful method, because translations
               do not make much sense out of the context of
               their transcript.  It may be better to use the
               TranscriptAdaptor::fetch_all_by_external_name instead.
               SQL wildcards % and _ are supported in the $external_name
               but their use is somewhat restricted for performance reasons.
               Users that really do want % and _ in the first three characters
               should use argument 3 to prevent optimisations
  Returntype : reference to a list of Translations
  Exceptions : none
  Caller     : general
  Status     : Medium Risk
             :   At some time may be deprecated to instead use 
             :   TranscriptAdaptor::fetch_all_by_external_name 

=cut
sub fetch_all_by_external_name {
  my ( $self, $external_name, $external_db_name, $override ) = @_;

  my $entry_adaptor = $self->db->get_DBEntryAdaptor();

  my @ids = $entry_adaptor->list_translation_ids_by_extids( 
            $external_name, $external_db_name, $override );

  my $transcript_adaptor = $self->db()->get_TranscriptAdaptor();

  my @out;
  foreach my $id (@ids) {
    my $transcript = $transcript_adaptor->fetch_by_translation_id($id);

    if ( defined($transcript) ) {
      push @out, $self->fetch_by_Transcript($transcript);
    }
  }

  return \@out;
}
*/

/* Don't bother
=head2 fetch_all_by_GOTerm

  Arg [1]   : Bio::EnsEMBL::OntologyTerm
              The GO term for which translations should be fetched.

  Example:  @translations = @{
              $translation_adaptor->fetch_all_by_GOTerm(
                $go_adaptor->fetch_by_accession('GO:0030326') ) };

  Description   : Retrieves a list of translations that are
                  associated with the given GO term, or with any of
                  its descendent GO terms.

  Return type   : listref of Bio::EnsEMBL::Translation
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
      $entryAdaptor->list_translation_ids_by_extids( $accession, 'GO' );
    foreach my $dbID (@ids) { $unique_dbIDs{$dbID} = 1 }
  }

  my @result;
  if ( scalar( keys(%unique_dbIDs) ) > 0 ) {
    my $transcript_adaptor = $self->db()->get_TranscriptAdaptor();

    foreach my $dbID ( sort { $a <=> $b } keys(%unique_dbIDs) ) {
      my $transcript =
        $transcript_adaptor->fetch_by_translation_id($dbID);
      if ( defined($transcript) ) {
        push( @result, $self->fetch_by_Transcript($transcript) );
      }
    }
  }

  return \@result;
} ## end sub fetch_all_by_GOTerm

=head2 fetch_all_by_GOTerm_accession

  Arg [1]   : String
              The GO term accession for which genes should be
              fetched.

  Example   :

    @genes =
      @{ $gene_adaptor->fetch_all_by_GOTerm_accession('GO:0030326') };

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

/* NIY
=head2 store

  Arg [1]    : Bio::EnsEMBL::Translation $translation
               The translation object to be stored in the database 
  Example    : $transl_id = $translation_adaptor->store($translation);
  Description: Stores a translation object in the database
  Returntype : int - the new dbID of the stored translation
  Exceptions : thrown if the dbID of the start_Exon or end_Exon is not 
               defined.
               thrown if only partial stable id information is present (e.g.
               identifier but not version number)
  Caller     : Transcript::store
  Status     : Stable

=cut
*/

IDType TranslationAdaptor_store(TranslationAdaptor *tlna, Translation *translation, IDType transcriptId) {
  
  Exon *startExon = Translation_getStartExon(translation);
  Exon *endExon   = Translation_getEndExon(translation);
 
  if (!startExon) {
    fprintf(stderr, "Translation must define a start_Exon to be stored.\n");
  }
 
  if (!endExon) {
    fprintf(stderr, "Translation must define an end_Exon to be stored.\n");
  }
 
  if (!Exon_getDbID(startExon)) {
    fprintf(stderr, "start_Exon must have a dbID for Translation to be stored.\n");
  }

  if (!Exon_getDbID(endExon)) {
    fprintf(stderr, "end_Exon must have a dbID for Translation to be stored.\n");
  }

  char qStr[1024];
  sprintf(qStr,
         "INSERT INTO translation " 
             "SET seq_start = %d,"
                " start_exon_id = "IDFMTSTR","
                " seq_end = %d," 
                " end_exon_id = "IDFMTSTR","
                " transcript_id = "IDFMTSTR,
          Translation_getStart(translation),
          Exon_getDbID(startExon),
          Translation_getEnd(translation),
          Exon_getDbID(endExon),
          transcriptId);

  if (Translation_getStableId(translation) != NULL) {
/*
      my $created = $self->db->dbc->from_seconds_to_date($translation->created_date());
      my $modified = $self->db->dbc->from_seconds_to_date($translation->modified_date());
*/
    // Assume version will be positive, Translation sets it to -1 when initialised
    int version = Translation_getVersion(translation) > 0 ? Translation_getVersion(translation) : 1; 
    sprintf(qStr,"%s, stable_id = '%s', version = %d, created_date = FROM_UNIXTIME(%ld), modified_date = FROM_UNIXTIME(%ld)",
            qStr, Translation_getStableId(translation), version, Translation_getCreated(translation), Translation_getModified(translation));
  }

  StatementHandle *sth = tlna->prepare((BaseAdaptor *)tlna,qStr,strlen(qStr));

  sth->execute(sth);
 
  IDType translDbID = sth->getInsertId(sth);

  sth->finish(sth);

  //
  // store object xref mappings to translations
  //
 
/* NIY
  DBEntryAdaptor *dbEntryAdaptor = DBAdaptor_getDBEntryAdaptor(tlna->dba);

  Vector *dbEntries = Translation_getAllDBEntries(translation);
  int i;
  for (i=0; i<Vector_getNumElement(dbEntries); i++) {
    DBEntry *dbe = Vector_getElementAt(dbEntries, i);
    DBEntryAdaptor_store(dbEntryAdaptor, dbe, translDbID, "Translation", 1);
  }
*/

  // storing the protein features associated with the translation
/* NIY
  my $pfadaptor = $self->db->get_ProteinFeatureAdaptor();
  foreach my $pf(@{$translation->get_all_ProteinFeatures}){
    $pfadaptor->store($pf, $transl_dbID);
  }
*/

// Slightly odd call in perl to get_all_Attributes when same call done in store call  $translation->get_all_Attributes();

  // store any translation attributes that are defined
  AttributeAdaptor *attrAdaptor = DBAdaptor_getAttributeAdaptor(tlna->dba);
  AttributeAdaptor_storeOnTranslationId(attrAdaptor, translDbID,
                                        Translation_getAllAttributes(translation, NULL));

  Translation_setDbID(translation, translDbID);
  Translation_setAdaptor(translation, (BaseAdaptor *)tlna);

  return translDbID;
}



/* NIY
=head2 remove

  Arg [1]    : Bio::EnsEMBL::Translation $translation
  Example    : $translation_adaptor->remove($translation);
  Description: Removes a translation completely from the database, and all
               associated information including protein features etc.
  Returntype : none
  Exceptions : throw on incorrect arguments
               warning if translation is not in this database
  Caller     : TranscriptAdaptor::remove
  Status     : Stable

=cut

sub remove {
  my $self = shift;
  my $translation = shift;

  if(!ref($translation) || !$translation->isa('Bio::EnsEMBL::Translation')) {
    throw("Bio::EnsEMBL::Translation argument expected.");
  }

  if( !$translation->is_stored($self->db()) ) {
    warning("Cannot remove translation " . $translation->dbID() . 
            ". Is not stored in this database.");
    return;
  }

  # remove athe attributes associated with this translation
  my $attrib_adp = $self->db->get_AttributeAdaptor;
  $attrib_adp->remove_from_Translation($translation);

  # remove all xref associations to this translation
  my $dbe_adaptor = $self->db()->get_DBEntryAdaptor();
  foreach my $dbe (@{$translation->get_all_DBEntries()}) {
    $dbe_adaptor->remove_from_object($dbe, $translation, 'Translation');
  }

  # remove all protein_features on this translation
  my $sth = $self->prepare
    ("DELETE FROM protein_feature WHERE translation_id = ?");
  $sth->bind_param(1,$translation->dbID,SQL_INTEGER);
  $sth->execute();
  $sth->finish();

  # remove the translation itself

  $sth = $self->prepare("DELETE FROM translation WHERE translation_id = ?" );
  $sth->bind_param(1,$translation->dbID,SQL_INTEGER);
  $sth->execute();
  $sth->finish();

  $translation->dbID( undef );
  $translation->adaptor(undef);

  return
}
*/


/*
=head2 list_dbIDs

  Arg [1]    : none
  Example    : @translation_ids = @{$translation_adaptor->list_dbIDs()};
  Description: Gets an array of internal ids for all translations in the current db
  Returntype : list of ints
  Exceptions : none
  Caller     : ?
  Status     : Stable

=cut
*/
Vector *TranslationAdaptor_listDbIDs(TranslationAdaptor *tlna, int ordered) {

  return TranslationAdaptor_listDbIDsPriv(tlna, "translation", NULL, ordered);
}

/*
=head2 list_stable_ids

  Arg [1]    : none
  Example    : @transl_stable_ids = @{$transl_adaptor->list_stable_dbIDs()};
  Description: Gets an array of stable ids for all translations in the current 
               db
  Returntype : reference to a list of strings
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/
Vector *TranslationAdaptor_listStableIDs(TranslationAdaptor *tlna) {

  return TranslationAdaptor_listDbIDsPriv(tlna, "translation", "stable_id", 0);
}

/*
=head2 _list_dbIDs

  Arg[1]      : String $table
  Arg[2]      : String $column
  Example     : $transl_adaptor->_list_dbIDs('translation','translation_id');
  Description : Local reimplementation to ensure multi-species translations
                are limited to their species alone
  Returntype  : ArrayRef of specified IDs
  Caller      : Internal
  Status      : Unstable
=cut
*/
// Hacky method 
Vector *TranslationAdaptor_listDbIDsPriv(TranslationAdaptor *tlna, char *table, char *column, int ordered) {
  Vector *ids;
  if (TranslationAdaptor_isMultiSpecies(tlna)) {
    fprintf(stderr,"multispecies db id fetch not implemented\n");
    exit(1);
/*
    $column ||= "${table}_id";
    my $sql = <<SQL;
select `tr`.`${column}` 
from translation tr
join transcript t using (transcript_id)
join seq_region sr using (seq_region_id)
join coord_system cs using (coord_system_id)
where cs.species_id =?
SQL
    return $self->dbc()->sql_helper()->execute_simple(-SQL => $sql, -PARAMS => [$self->species_id()]);
*/
  } else {
// NIY: Shouldn't really do this direct BaseAdaptor call but ...
    ids = BaseAdaptor_listDbIDs((BaseAdaptor *)tlna, table, column, ordered);
  }
  return ids;
}


/* Lets try to get away without this - it doesn't seem very necessary
=head2 fetch_by_dbID

  Arg [1]    : int $dbID
               The internal identifier of the Translation to obtain
  Example    : $translation = $translation_adaptor->fetch_by_dbID(1234);
  Description: This fetches a Translation object via its internal id.
               This is only debatably useful since translations do
               not make much sense outside of the context of their
               Transcript.  Consider using fetch_by_Transcript instead.
  Returntype : Bio::EnsEMBL::Translation, or undef if the translation is not
               found.
  Exceptions : warning if an additional (old style) Transcript argument is
               provided
  Caller     : ?
  Status     : Stable

=cut
sub fetch_by_dbID {
  my ( $self, $dbID, $transcript ) = @_;

  if ($transcript) {
    deprecate(   "Use of fetch_by_dbID "
               . "with a Transcript argument is deprecated."
               . "Use fetch_by_Transcript instead." );
  }

  if ( !defined($dbID) ) {
    throw("dbID argument is required");
  }

  my $transcript_adaptor = $self->db()->get_TranscriptAdaptor();
  $transcript = $transcript_adaptor->fetch_by_translation_id($dbID);

  if ( defined($transcript) ) {
    my $translation = $self->fetch_by_Transcript($transcript);

    if ( defined($translation) && $translation->dbID()==$dbID ) {
      return $translation;
    }

    my @alt_translations =
      @{ $self->fetch_all_alternative_by_Transcript($transcript) };

    foreach my $alt_translation (@alt_translations) {
      if ( $alt_translation->dbID() == $dbID ) {
        return $alt_translation;
      }
    }
  }

  return undef;
}
*/

/*
=head2 fetch_by_stable_id

  Arg [1]    : string $stable_id
               The stable identifier of the Translation to obtain
  Example    : $translation = $translation_adaptor->fetch_by_stable_id("ENSP00001");
  Description: This fetches a Translation object via its stable id.
               This is only debatably useful since translations do
               not make much sense outside of the context of their
               Transcript.  Consider using fetch_by_Transcript instead.
  Returntype : Bio::EnsEMBL::Translation or undef if the translation is not
               found.
  Exceptions : warning if an additional (old style) Transcript argument is
               provided
  Caller     : ?
  Status     : Stable

=cut
*/
Translation *TranslationAdaptor_fetchByStableId(TranslationAdaptor *tlna, char *stableId) {
  if (stableId == NULL) {
    fprintf(stderr, "stable id argument is required\n");
    exit(1);
  }

  TranscriptAdaptor *ta = DBAdaptor_getTranscriptAdaptor(tlna->dba);

  Transcript *transcript = TranscriptAdaptor_fetchByTranslationStableId(ta, stableId);

  if (!transcript) {
    return NULL;
  } else {
    return TranslationAdaptor_fetchByTranscript(tlna, transcript);
  }
}


/*
=head2 fetch_all_by_Transcript_list

  Arg [1]    : reference to list of Bio::EnsEMBL::Transcripts $transcripts
               The list of $transcripts to obtain Translation object for.
  Example    : @translations = @{$tla->fetch_all_by_Transcript_list([$t1,$t2]);
  Description: Fetches all translations associated with the list of transcripts
               passed to this method.  The passed transcripts will also have
               their translation set by this method.
  Returntype : Reference to list of Bio::EnsEMBL::Translations
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut
*/
Vector *TranslationAdaptor_fetchAllByTranscriptList(TranslationAdaptor *tlna, Vector *transcripts) {
  if (transcripts == NULL) {
    fprintf(stderr,"Need array of Transcripts as argument in TranslationAdaptor_fetchAllByTranscriptList - bye");
    exit(1);
  }

  if (Vector_getNumElement(transcripts) == 0) {
    return Vector_new();
  }

  IDHash *transHash = IDHash_new(IDHASH_MEDIUM);
  int i;
  for (i=0; i<Vector_getNumElement(transcripts); i++) {
    Transcript *trans = Vector_getElementAt(transcripts, i);
    IDType id = Transcript_getDbID(trans);

    if (!IDHash_contains(transHash, id)) {
      IDHash_add(transHash, id, trans);
    }
  }
  IDType *uniqueIds = IDHash_getKeys(transHash);
  int nUniqueId = IDHash_getNumValues(transHash);

  Vector *out = Vector_new();

  // mysql is faster ?? really ?? and we ensure that we do not exceed the max query size by
  // splitting large queries into smaller queries of 200 ids
  // Note: The exon fetching code in exon adaptor doesn't do this splitting, so any attempted
  // maximum length checking here is basically useless
  int maxSize = 200;

  // Unused in perlmy %ex_hash;

  for (i=0; i<nUniqueId; i+=maxSize) {
    char idStr[655500];
    idStr[0] = '\0';

    // Special case for one remaining Id
    if (i == nUniqueId-1) {
      sprintf(idStr, " = "IDFMTSTR, uniqueIds[i]);
    } else {
      char tmpStr[1024];
      strcat(idStr, " IN (");
      int j;
      for (j=0; j<maxSize && j+i<nUniqueId; j++) {
        if (j!=0) {
          strcat(idStr,", ");
        }
        sprintf(tmpStr, IDFMTSTR, uniqueIds[i+j]);
        strcat(idStr, tmpStr);
      }
      strcat(idStr, ")");
    }
    
    IDHash *canonicalLookup = IDHash_new(IDHASH_SMALL);
    // Annoying new style code - have to guess a bit what it does
    //my $canonical_lookup = $self->dbc()->sql_helper()->execute_into_hash(
    //  -SQL => 'SELECT transcript_id, canonical_translation_id FROM transcript WHERE transcript_id '.$id_str
    //);
    
    char qStr[655500];
    sprintf(qStr,"SELECT transcript_id, canonical_translation_id FROM transcript WHERE transcript_id %s", idStr);
    StatementHandle *sth = tlna->prepare((BaseAdaptor *)tlna,qStr,strlen(qStr));
    sth->execute(sth);
    ResultRow *row;
    while ((row = sth->fetchRow(sth))) {
      IDType transcriptId           = row->getLongLongAt(row,0);
      IDType *idP; 
      IDType canonicalTranslationId = row->getLongLongAt(row,1);

      if (!IDHash_contains(canonicalLookup, transcriptId)) {
        if ((idP = calloc(1,sizeof(IDType))) == NULL) {
          fprintf(stderr, "Failed allocating space for a id\n");
          exit(1);
        }

        *idP = canonicalTranslationId;

        IDHash_add(canonicalLookup, transcriptId, idP);
      }
    }
    sth->finish(sth);

    // Just use UNIX_TIMESTAMP
    // Can't be arsed my $created_date = $self->db->dbc->from_date_to_seconds("tl.created_date");
    // Can't be arsed my $modified_date = $self->db->dbc->from_date_to_seconds("tl.modified_date");

    sprintf(qStr,"SELECT tl.transcript_id, tl.translation_id, tl.start_exon_id, "
                        "tl.end_exon_id, tl.seq_start, tl.seq_end, "
                        "tl.stable_id, tl.version, UNIX_TIMESTAMP(tl.created_date), UNIX_TIMESTAMP(tl.modified_date) "
                   "FROM translation tl "
                  "WHERE tl.transcript_id %s", idStr);

    sth = tlna->prepare((BaseAdaptor *)tlna,qStr,strlen(qStr));
    sth->execute(sth);

    while ((row = sth->fetchRow(sth))) {
      // Need transcriptId so fetch it here as well as in translationFromResultRow
      IDType transcriptId  = row->getLongLongAt(row,0); 
      Transcript *tr = IDHash_getValue(transHash, transcriptId);

      Translation *tl = TranslationAdaptor_translationFromResultRow(tlna, row, tr);
  
      IDType canonicalTranslationId = *((IDType *)IDHash_getValue(canonicalLookup, Transcript_getDbID(tr)));
      if (Translation_getDbID(tl) == canonicalTranslationId) {
        Transcript_setTranslation(tr, tl);
      }

      Vector_addElement(out, tl);
    }
    sth->finish(sth);
    IDHash_free(canonicalLookup, free);
  }

  free(uniqueIds);
  IDHash_free(transHash,NULL);

  return out;
}

// This avoids code duplication in several methods, when fetching translation data from the db
Translation *TranslationAdaptor_translationFromResultRow(TranslationAdaptor *tlna, ResultRow *row, Transcript *transcript) {

  IDType transcriptId  = row->getLongLongAt(row,0); 
  IDType translationId = row->getLongLongAt(row,1);
  IDType startExonId   = row->getLongLongAt(row,2);
  IDType endExonId     = row->getLongLongAt(row,3);
  long seqStart        = row->getLongAt(row,4);
  long seqEnd          = row->getLongAt(row,5);
  char *stableId       = row->getStringAt(row,6);
  int version          = row->getIntAt(row,7);
  int createdDate      = row->getIntAt(row,8);
  int modifiedDate     = row->getIntAt(row,9);

  
  // TranslationAdaptor_fetchAllAlternativeByTranscript had test for undef translationId. Equivalent will be 0 in C (ints can't be undef, and 0 translationId will be invalid in Db)
  
  if (translationId == 0) {
    return NULL;
  }

  Exon *startExon = NULL;
  Exon *endExon = NULL;

  // this will load all the exons whenever we load the translation
  // but I guess thats ok ....

  int i;
  for (i=0; i<Transcript_getExonCount(transcript); i++) {
    Exon *exon = Transcript_getExonAt(transcript,i);

    if(startExon == NULL && Exon_getDbID(exon) == startExonId ) {
      startExon = exon;
      if (endExon != NULL) break;
    }

    if(endExon == NULL && Exon_getDbID(exon) == endExonId ) {
      endExon = exon;
      if (startExon != NULL) break;
    }
  }

  if (startExon == NULL || endExon == NULL) {
    fprintf(stderr,"Could not find start or end exon in transcript\n");
    exit(1);
  }

  Translation *tl = Translation_new();
  Translation_setDbID(tl, translationId);
  Translation_setStart(tl, seqStart);
  Translation_setEnd(tl, seqEnd);
  Translation_setStartExon(tl, startExon);
  Translation_setEndExon(tl, endExon);
  Translation_setStableId(tl, stableId);
  Translation_setVersion(tl, version);
// Note dates had undef option in perl -created_date => $created_date || undef,
  Translation_setCreated(tl, createdDate);
  Translation_setModified(tl, modifiedDate);
  Translation_setAdaptor(tl, (BaseAdaptor *)tlna);

  return tl;
}


/*
=head2 fetch_all_by_DBEntry

  Description: DEPRECATED, this has been renames fetch_all_by_external_name

=cut

sub fetch_all_by_DBEntry {
  my $self = shift;
  deprecate("Use fetch_all_by_external_name instead.");
  return $self->fetch_all_by_external_name(@_);
}

=head2 get_stable_entry_info

 Description: DEPRECATED - This method should no longer be needed. Stable
              id info is fetched when the transcript is.

=cut

sub get_stable_entry_info {
  my ($self,$translation) = @_;

  deprecate( "This method shouldnt be necessary any more" );

  unless(defined $translation && ref $translation && 
	 $translation->isa('Bio::EnsEMBL::Translation') ) {
    throw("Needs a Translation object, not a [$translation]");
  }

  my $sth = $self->prepare("SELECT stable_id, version 
                            FROM   translation
                            WHERE  translation_id = ?");
  $sth->bind_param(1,$translation->dbID,SQL_INTEGER);
  $sth->execute();

  my @array = $sth->fetchrow_array();
  $translation->{'_stable_id'} = $array[0];
  $translation->{'_version'}   = $array[1];

  return 1;
}
*/

/*
=head2 fetch_all

  Example     : $translations = $translation_adaptor->fetch_all();
  Description : Retrieves all canonical and alternative translations 
                stored in the database.
  Returntype  : listref of Bio::EnsEMBL::Translation
  Caller      : general
  Status      : At Risk

=cut
*/
Vector *TranslationAdaptor_fetchAll(TranslationAdaptor *tlna) {
  TranscriptAdaptor *ta = DBAdaptor_getTranscriptAdaptor(tlna->dba);

  Vector *translations;
  Vector *transcripts = TranscriptAdaptor_fetchAll(ta);

  int i;
  for (i=0; i<Vector_getNumElement(transcripts); i++) {
    Transcript *transcript = Vector_getElementAt(transcripts, i);
  
    Translation *translation = TranslationAdaptor_fetchByTranscript(tlna, transcript);
    if (translation != NULL) {
      Vector_addElement(translations, translation);
    }
    Vector *altTranslations = TranslationAdaptor_fetchAllAlternativeByTranscript(tlna, transcript);
    int j;
    for (j=0; j<Vector_getNumElement(altTranslations); j++) {
      Vector_addElement(translations, Vector_getElementAt(altTranslations, j));
    }
    Vector_free(altTranslations);
  }
// NIY: Do I need to free transcripts properly??
  Vector_free(transcripts);
  return translations;
}

