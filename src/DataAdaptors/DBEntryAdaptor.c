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

#include "DBEntryAdaptor.h"
#include "DBAdaptor.h"
#include "DBEntry.h"
#include "IDHash.h"
#include "Class.h"
#include "RawContig.h"
#include "Transcript.h"
#include "Gene.h"
#include "Translation.h"

Vector *DBEntryAdaptor_fetchByObjectType(DBEntryAdaptor *dbea, IDType ensObj, char *ensType);


DBEntryAdaptor *DBEntryAdaptor_new(DBAdaptor *dba) {
  DBEntryAdaptor *dea;

  if ((dea = (DBEntryAdaptor *)calloc(1,sizeof(DBEntryAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for DBEntryAdaptor\n");
    return NULL;
  }
  BaseAdaptor_init((BaseAdaptor *)dea, dba, DBENTRY_ADAPTOR);

  return dea;
}

DBEntry *DBEntryAdaptor_fetchByDbID(DBEntryAdaptor *dbea, IDType dbID) {
  char qStr[512];
  StatementHandle *sth;
  ResultRow *row;
  DBEntry *dbe = NULL;

  sprintf(qStr,
   "SELECT xref.xref_id, xref.dbprimary_acc, xref.display_label,"
   "       xref.version, xref.description,"
   "       exDB.db_name, exDB.db_release, es.synonym"
   " FROM  (xref, external_db exDB)"
   " LEFT JOIN external_synonym es on es.xref_id = xref.xref_id"
   " WHERE  xref.xref_id = " IDFMTSTR
   " AND    xref.external_db_id = exDB.external_db_id",
   dbID);
   

  sth = dbea->prepare((BaseAdaptor *)dbea,qStr,strlen(qStr));

  sth->execute(sth);
  

  // Why???? my %duplicate;

  while ((row = sth->fetchRow(sth))){
    if (!row->col(row,0)) {
      fprintf(stderr,"WARNING: Got xref with no refID\n");
      return NULL;
    }
    
    if (!dbe) {
      dbe = DBEntry_new();

      DBEntry_setAdaptor(dbe,(BaseAdaptor *)dbea);
      DBEntry_setDbID(dbe, dbID);
      DBEntry_setPrimaryId(dbe, row->getStringAt(row,1));
      DBEntry_setDisplayId(dbe, row->getStringAt(row,2));
      DBEntry_setVersion(dbe, row->getStringAt(row,3));
      DBEntry_setRelease(dbe, row->getStringAt(row,6));
      DBEntry_setDbName(dbe, row->getStringAt(row,5));
      
      if (row->col(row,4)) DBEntry_setDescription(dbe, row->getStringAt(row,4));
    }
 
    if (row->col(row,7)) DBEntry_addSynonym(dbe, row->getStringAt(row,7));
  }

  sth->finish(sth);

  return dbe;
}


IDType DBEntryAdaptor_store(DBEntryAdaptor *dbea, DBEntry *exObj, 
                            IDType ensObject, char *ensType, int ignoreRelease) {
  fprintf(stderr,"DBEntryAdaptor_store does not implement ignoreRelease functionality yet\n");

  char qStr[512];
  StatementHandle *sth;
  ResultRow *row;
  IDType dbRef;
  IDType dbX;

  //
  // Check for the existance of the external_db, throw if it does not exist
  //
  sprintf(qStr,
     "SELECT external_db_id"
     "  FROM external_db"
     " WHERE db_name = '%s'"
     "   AND db_release = %s",
     DBEntry_getDbName(exObj),
     DBEntry_getRelease(exObj));

  sth = dbea->prepare((BaseAdaptor *)dbea,qStr,strlen(qStr));
  sth->execute(sth);
    
  row = sth->fetchRow(sth);
  if( row == NULL ) {
    sth->finish(sth);
    fprintf(stderr,"Error: external_db [%s] release [%s] does not exist\n", 
            DBEntry_getDbName(exObj), DBEntry_getRelease(exObj));
    exit(1);
  }

  dbRef =  row->getLongLongAt(row,0);
  sth->finish(sth);
    
  //
  // Check for the existance of the external reference, add it if not present
  //
  sprintf(qStr,
       "SELECT xref_id"
       "  FROM xref"
       " WHERE external_db_id = " IDFMTSTR
       "   AND dbprimary_acc = '%s'"
       "   AND version = %s",
      dbRef,
      DBEntry_getPrimaryId(exObj),
      DBEntry_getVersion(exObj));

  sth = dbea->prepare((BaseAdaptor *)dbea,qStr,strlen(qStr));
  sth->execute(sth);

  row = sth->fetchRow(sth);
    
  if (row != NULL) {
    dbX =  row->getLongLongAt(row,0);
    sth->finish(sth);
  } else {
    //
    // store the new xref
    //

    // First finish the old sth
    sth->finish(sth);

// NIY Handling NULL values
    sprintf(qStr,
       "INSERT ignore INTO xref"
       " SET dbprimary_acc = '%s',"
       "    display_label = '%s',"
       "    version = %s,"
       "    description = '%s',"
       "    external_db_id = " IDFMTSTR,
       DBEntry_getPrimaryId(exObj),
       DBEntry_getDisplayId(exObj),
       DBEntry_getVersion(exObj),
       DBEntry_getDescription(exObj),
       dbRef
      );
    sth = dbea->prepare((BaseAdaptor *)dbea,qStr,strlen(qStr));

    sth->execute(sth);
    dbX = sth->getInsertId(sth);

    sth->finish(sth);
	
    //
    // store the synonyms for the new xref
    // 
    if (DBEntry_getAllSynonyms(exObj)) {
      StatementHandle *checkSth;
      StatementHandle *storeSth;
      int i;
      Vector *synonyms;

      sprintf(qStr,
              "SELECT xref_id, synonym"
              " FROM external_synonym"
              " WHERE xref_id = %" IDFMTSTR
              " AND synonym = '%%s'");

      checkSth = dbea->prepare((BaseAdaptor *)dbea,qStr,strlen(qStr));

      sprintf(qStr,
        "INSERT ignore INTO external_synonym"
        " SET xref_id = %" IDFMTSTR ", synonym = '%%s'");     

      storeSth = dbea->prepare((BaseAdaptor *)dbea,qStr,strlen(qStr));

      synonyms = DBEntry_getAllSynonyms(exObj);

      for (i=0;i<Vector_getNumElement(synonyms); i++) {	    
        char *syn = Vector_getElementAt(synonyms,i);
        checkSth->execute(checkSth, dbX, syn);
        row = checkSth->fetchRow(checkSth);
        if (!row) {
          storeSth->execute(storeSth, dbX, syn);
        }
      }
  	
      checkSth->finish(checkSth);
      storeSth->finish(storeSth);
    }
  }

  //
  // check if the object mapping was already stored
  //
  sprintf(qStr,
           "SELECT xref_id"
           " FROM object_xref"
           " WHERE xref_id = " IDFMTSTR
           " AND   ensembl_object_type = '%s'"
           " AND   ensembl_id = " IDFMTSTR,
         dbX, ensType, ensObject);

  sth = dbea->prepare((BaseAdaptor *)dbea,qStr,strlen(qStr));

  sth->execute(sth);

  row = sth->fetchRow(sth);
// NOTE row will be invalid after this call but will still
//      indicate whether something was found
  sth->finish(sth);
    
  if (!row) {
    IDType Xidt;

    //
    // Store the reference to the internal ensembl object
    //
    sprintf(qStr,
         "INSERT ignore INTO object_xref"
         " SET xref_id = " IDFMTSTR ","
         "     ensembl_object_type = '%s',"
         "     ensembl_id = " IDFMTSTR,
        dbX, ensType, ensObject);

    sth = dbea->prepare((BaseAdaptor *)dbea,qStr,strlen(qStr));
	
    sth->execute(sth);
    DBEntry_setDbID(exObj, dbX);
    DBEntry_setAdaptor(exObj, (BaseAdaptor *)dbea);
      
    Xidt = sth->getInsertId(sth);

    //
    // If this is an IdentityXref need to store in that table too
    //
    if (DBEntry_getIdentityXref(exObj)) {
      IdentityXref *idx = DBEntry_getIdentityXref(exObj);
      sprintf(qStr,
             "INSERT ignore INTO identity_xref"
             " SET object_xref_id = " IDFMTSTR ","
             "     query_identity = %f,"
             "     target_identity = %f",
             Xidt, 
             IdentityXref_getQueryIdentity(idx),
             IdentityXref_getTargetIdentity(idx));

      sth = dbea->prepare((BaseAdaptor *)dbea,qStr,strlen(qStr));
      sth->execute(sth);
      sth->finish(sth);
    }
  } 
  return dbX;    
}

IDType DBEntryAdaptor_exists(DBEntryAdaptor *dbea, DBEntry *dbe) {
  char qStr[512];
  StatementHandle *sth;
  ResultRow *row;
  IDType dbID;

  if (!dbe || !Class_isDescendent(CLASS_DBENTRY, dbe->objectType)) {
    fprintf(stderr,"Error: arg must be a DBEntry\n");
    exit(1);
  }

// NIY Was dbe->external_db instead of dbe->dbname - are they different?
// NIY mysql_quote strings
  sprintf(qStr,
              "SELECT x.xref_id "
              " FROM   xref x, external_db xdb"
              " WHERE  x.external_db_id = xdb.external_db_id"
              " AND    x.display_label = '%s' "
              " AND    xdb.db_name = '%s'",
          DBEntry_getDisplayId(dbe),
          DBEntry_getDbName(dbe));
  

  sth = dbea->prepare((BaseAdaptor *)dbea,qStr,strlen(qStr));

  sth->execute(sth);

  row = sth->fetchRow(sth);
  if( row == NULL ) {
    sth->finish(sth);
    return 0;
  }

  dbID = row->getLongLongAt(row,0);

  sth->finish(sth);

  return dbID;
}

int DBEntryAdaptor_fetchAllByGene(DBEntryAdaptor *dbea, Gene *gene) {
  char qStr[512];
  StatementHandle *sth;
  ResultRow *row;
 
  sprintf(qStr,
      "SELECT t.transcript_id, t.canonical_translation_id"
      " FROM   transcript t"
      " WHERE  t.gene_id = " IDFMTSTR,
       Gene_getDbID(gene));

  sth = dbea->prepare((BaseAdaptor *)dbea,qStr,strlen(qStr));
  
  sth->execute(sth);
  
  while ((row = sth->fetchRow(sth))) {
    IDType transcriptId = row->getLongLongAt(row,0);
    int i;
    Vector *transLinks;

    if (row->col(row,1)) {
      IDType translationId = row->getLongLongAt(row,1);
      Vector *translatLinks = DBEntryAdaptor_fetchByObjectType(dbea, translationId,"Translation");

      for (i=0;i<Vector_getNumElement(translatLinks); i++) {
        Gene_addDBLink(gene,Vector_getElementAt(translatLinks,i));
      }
      Vector_free(translatLinks);
    }

    
    transLinks = DBEntryAdaptor_fetchByObjectType(dbea, transcriptId,"Transcript");
    for (i=0;i<Vector_getNumElement(transLinks); i++) {
      Gene_addDBLink(gene, Vector_getElementAt(transLinks,i));
    }
    Vector_free(transLinks);
  }

/* NIY This is wrong so I'm not going to implement it!
  if($gene->stable_id){
    my $genelinks = $self->_fetch_by_object_type( $gene->stable_id, 'Gene' );
    foreach my $genelink ( @$genelinks ) {
      $gene->add_DBLink( $genelink );
    }
  }
*/
  return 1;
}

Vector *DBEntryAdaptor_fetchAllByRawContig(DBEntryAdaptor *dbea, RawContig *contig) {
  Vector *contigLinks = DBEntryAdaptor_fetchByObjectType(dbea, RawContig_getDbID(contig),"RawContig");
  return contigLinks;
}

int DBEntryAdaptor_fetchAllByTranscript(DBEntryAdaptor *dbea, Transcript *trans) {
  char qStr[512];
  StatementHandle *sth;
  ResultRow *row;
  Vector *transLinks;
  int i;

  sprintf(qStr, 
    "SELECT t.canonical_translation_id" 
    " FROM transcript t"
    " WHERE t.transcript_id = " IDFMTSTR,
    Transcript_getDbID(trans));
  
  sth = dbea->prepare((BaseAdaptor *)dbea,qStr,strlen(qStr));

  sth->execute(sth);

  // 
  // Did this to be consistent with fetch_by_Gene, but don't like
  // it (filling in the object). I think returning the array would
  // be better. Oh well. EB
  //
  
  while ((row = sth->fetchRow(sth))) {
    IDType translationId = row->getLongLongAt(row,0);
    Vector *translatLinks = DBEntryAdaptor_fetchByObjectType(dbea, translationId,"Translation");
    for (i=0;i<Vector_getNumElement(translatLinks); i++) {
      Transcript_addDBLink(trans,Vector_getElementAt(translatLinks,i));
    }
    Vector_free(translatLinks);
  }

  sth->finish(sth);

  transLinks = DBEntryAdaptor_fetchByObjectType(dbea, Transcript_getDbID(trans),"Transcript");
      fprintf(stderr,"transLinks\n");
  for (i=0;i<Vector_getNumElement(transLinks); i++) {
    Transcript_addDBLink(trans,Vector_getElementAt(transLinks,i));
  }
  Vector_free(transLinks);

  return 1;
}

Vector *DBEntryAdaptor_fetchAllByTranslation(DBEntryAdaptor *dbea, Translation *trans) {
  Vector *translatLinks = DBEntryAdaptor_fetchByObjectType(dbea, Translation_getDbID(trans),"Translation");
  return translatLinks;
}

Vector *DBEntryAdaptor_fetchByObjectType(DBEntryAdaptor *dbea, IDType ensObj, char *ensType) {
  Vector *out;
  char qStr[1024];
  StatementHandle *sth;
  ResultRow *row;
  IDHash *seen;
  
  if (!ensObj) {
    fprintf(stderr,"Error: Can't fetchByObjectType without an object\n");
    exit(1);
  }

  if (!ensType) {
    fprintf(stderr,"Error: Can't fetchByObjectType without a type\n");
    exit(1);
  }

// Not sure if idt identities are right way round
  sprintf(qStr,
    "SELECT xref.xref_id, xref.dbprimary_acc, xref.display_label, xref.version,"
    "       xref.description,"
    "       exDB.db_name, exDB.db_release, exDB.status," 
    "       oxr.object_xref_id,"
    "       es.synonym," 
    "       idt.xref_identity, idt.ensembl_identity"
    " FROM  (external_db exDB, object_xref oxr, xref xref)" 
    " LEFT JOIN external_synonym es on es.xref_id = xref.xref_id"
    " LEFT JOIN identity_xref idt on idt.object_xref_id = oxr.object_xref_id"
    " WHERE  xref.xref_id = oxr.xref_id"
    "  AND  xref.external_db_id = exDB.external_db_id"
    "  AND  oxr.ensembl_id = " IDFMTSTR
    "  AND  oxr.ensembl_object_type = '%s'",
    ensObj,
    ensType);
  
  sth = dbea->prepare((BaseAdaptor *)dbea,qStr,strlen(qStr));

  sth->execute(sth);

  seen = IDHash_new(IDHASH_SMALL);
  out = Vector_new();

  while ((row = sth->fetchRow(sth))) {
    DBEntry *exDB;
    IDType refID = row->getLongLongAt(row,0);
			    
    // using an outer join on the synonyms as well as on identity_xref, we
    // now have to filter out the duplicates (see v.1.18 for
    // original). Since there is at most one identity_xref row per xref,
    // this is easy enough; all the 'extra' bits are synonyms

    if (!IDHash_contains(seen,refID))  {
      exDB = DBEntry_new();
      DBEntry_setAdaptor(exDB,(BaseAdaptor *)dbea);
      DBEntry_setDbID(exDB, refID);
      DBEntry_setPrimaryId(exDB, row->getStringAt(row,1));
      DBEntry_setDisplayId(exDB, row->getStringAt(row,2));
      DBEntry_setVersion(exDB, row->getStringAt(row,3));
      DBEntry_setDbName(exDB, row->getStringAt(row,5));
      DBEntry_setRelease(exDB, row->getStringAt(row,6));

      if (row->col(row,10)) {
        IdentityXref *idx = IdentityXref_new();
        DBEntry_setIdentityXref(exDB,idx);
	IdentityXref_setQueryIdentity(idx, row->getDoubleAt(row,10));
	IdentityXref_setTargetIdentity(idx, row->getDoubleAt(row,11));
      }
      
      if (row->col(row,4)) DBEntry_setDescription(exDB, row->getStringAt(row,4));
      if (row->col(row,7)) DBEntry_setStatus(exDB, row->getStringAt(row,7));
      
      Vector_addElement(out, exDB);
      IDHash_add(seen, refID, exDB);
    } 

    exDB = IDHash_getValue(seen, refID);

    if (row->col(row,9)) {
      DBEntry_addSynonym(exDB,row->getStringAt(row,9));
    }
  }

  IDHash_free(seen, NULL);

  sth->finish(sth);
  
  return out;
}

#ifdef DONE
/*
=head2 list_gene_ids_by_extids

  Arg [1]    : string $external_id
  Example    : none
  Description: Retrieve a list of geneid by an external identifier that is linked to 
               any of the genes transcripts, translations or the gene itself 
  Returntype : listref of strings
  Exceptions : none
  Caller     : unknown

=cut

sub list_gene_ids_by_extids{
   my ($self,$name) = @_;

   my %T = map { ($_,1) }
       $self->_type_by_external_id( $name, 'Translation', 'gene' ),
       $self->_type_by_external_id( $name, 'Transcript',  'gene' ),
       $self->_type_by_external_id( $name, 'Gene' );
   return keys %T;
}
=head2 list_transcript_ids_by_extids

  Arg [1]    : string $external_id
  Example    : none
  Description: Retrieve a list transcriptid by an external identifier that is linked to 
               any of the genes transcripts, translations or the gene itself 
  Returntype : listref of strings
  Exceptions : none
  Caller     : unknown

=cut

sub list_transcript_ids_by_extids{
   my ($self,$name) = @_;
   my @transcripts;

   my %T = map { ($_,1) }
       $self->_type_by_external_id( $name, 'Translation', 'transcript' ),
       $self->_type_by_external_id( $name, 'Transcript' );
   return keys %T;
}

=head2 list_translation_ids_by_extids

  Arg [1]    :  string $name 
  Example    :  none
  Description:  Gets a list of translation IDs by external display IDs
  Returntype :  list of Ints
  Exceptions :  none
  Caller     :  unknown

=cut

sub list_translation_ids_by_extids{
  my ($self,$name) = @_;
  return $self->_type_by_external_id( $name, 'Translation' );
}

=head2 _type_by_external_id

  Arg [1]    : string $name
  			   (dbprimary_acc)
  Arg [2]    : string $ensType
  			   (Object_type)
  Arg [3]    : string $extraType
  			   (other object type to be returned) - optional
  Example    : $self->_type_by_external_id( $name, 'Translation' ) 
  Description: Gets
  Returntype : list of ensembl_IDs
  Exceptions : none
  Caller     : list_translation_ids_by_extids
               translationids_by_extids
  			   geneids_by_extids

=cut


sub _type_by_external_id{
  my ($self,$name,$ensType,$extraType) = @_;
   
  my $from_sql = '';
  my $where_sql = '';
  my $ID_sql = "oxr.ensembl_id";
  if(defined $extraType) {
    $ID_sql = "t.${extraType}_id";
    $from_sql = 'transcript as t, ';
    $where_sql = 't.'.lc($ensType).'_id = oxr.ensembl_id and ';
  }
  my @queries = (
    "select $ID_sql
       from $from_sql xref, object_xref as oxr
      where $where_sql xref.dbprimary_acc = ? and
            xref.xref_id = oxr.xref_id and oxr.ensembl_object_type= ?",
    "select $ID_sql 
       from $from_sql xref, object_xref as oxr
      where $where_sql xref.display_label = ? and
            xref.xref_id = oxr.xref_id and oxr.ensembl_object_type= ?",
    "select $ID_sql
       from $from_sql object_xref as oxr, external_synonym as syn
      where $where_sql syn.synonym = ? and
            syn.xref_id = oxr.xref_id and oxr.ensembl_object_type= ?",
  );

// Increase speed of query by splitting the OR in query into three separate 
// queries. This is because the 'or' statments render the index useless 
// because MySQL can't use any fields in the index.
  
  my %hash = (); 
  foreach( @queries ) {
    my $sth = $self->prepare( $_ );
  sth = dbea->prepare((BaseAdaptor *)dbea,qStr,strlen(qStr));
    $sth->execute("$name", $ensType);
    while( my $r = $sth->fetchrow_array() ) {
      $hash{$r} = 1;
    }
  }
  return keys %hash;
}
*/
#endif


