#include "TranscriptAdaptor.h"
#include "ExonAdaptor.h"
#include "BaseAdaptor.h"
#include "MysqlUtil.h"
#include "TranslationAdaptor.h"

#include "StatementHandle.h"
#include "ResultRow.h"

TranscriptAdaptor *TranscriptAdaptor_new(DBAdaptor *dba) {
  TranscriptAdaptor *ta;

  if ((ta = (TranscriptAdaptor *)calloc(1,sizeof(TranscriptAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for TranscriptAdaptor\n");
    return NULL;
  }
  BaseAdaptor_init((BaseAdaptor *)ta, dba, TRANSCRIPT_ADAPTOR);

  return ta;
}

Transcript *TranscriptAdaptor_fetchByDbID(TranscriptAdaptor *ta, IDType dbID) {
  Transcript *trans;
  char qStr[256];
  StatementHandle *sth;
  ResultRow *row;
  ExonAdaptor *ea;
  int seen = 0;

  trans = Transcript_new();
  ea = DBAdaptor_getExonAdaptor(ta->dba);

  sprintf(qStr,
    "SELECT exon_id"
    " FROM   exon_transcript"
    " WHERE  transcript_id = "
    IDFMTSTR
    " ORDER BY rank",dbID);

  sth = ta->prepare((BaseAdaptor *)ta,qStr,strlen(qStr));
  sth->execute(sth);

  while ((row = sth->fetchRow(sth))) {
    Exon *exon = ExonAdaptor_fetchByDbID(ea, row->getLongLongAt(row,0));
    Transcript_addExon(trans,exon);
    seen = 1;
  }
  sth->finish(sth);

  if (!seen) {
    fprintf(stderr,"ERROR: transcript " IDFMTSTR "is not present in db",dbID);
    return NULL;
  }
  Transcript_setDbID(trans, dbID);
  Transcript_setAdaptor(trans, (BaseAdaptor *)ta );

  sprintf(qStr,
    "SELECT translation_id"
    " FROM  transcript"
    " WHERE  transcript_id = " IDFMTSTR,dbID);

  sth = ta->prepare((BaseAdaptor *)ta,qStr,strlen(qStr));
  sth->execute(sth);

  row = sth->fetchRow(sth);
  if( row != NULL ) {
    Transcript_setTranslationId(trans, row->getLongLongAt(row,0));
  }
  sth->finish(sth);

  return trans;
}

int TranscriptAdaptor_getStableEntryInfo(TranscriptAdaptor *ta, Transcript *transcript) {
  char qStr[256];
  StatementHandle *sth;
  ResultRow *row;

  if( !transcript ) {
    fprintf(stderr, "ERROR: TranscriptAdaptor_getStableEntryInfo needs a transcript object\n");
    exit(1);
  }

  sprintf(qStr,
          "SELECT stable_id, version"
          " FROM transcript_stable_id"
          " WHERE transcript_id = " IDFMTSTR, Transcript_getDbID(transcript));

  sth = ta->prepare((BaseAdaptor *)ta,qStr,strlen(qStr));
  sth->execute(sth);

  row = sth->fetchRow(sth);
  if( row == NULL ) {
    fprintf(stderr,"WARNING: Failed fetching stable id info\n");
    sth->finish(sth);
    return 0;
  }

  Transcript_setStableId(transcript,row->getStringAt(row,0));
  Transcript_setVersion(transcript,row->getIntAt(row,1));

  sth->finish(sth);

  return 1;
}

IDType TranscriptAdaptor_store(TranscriptAdaptor *ta, Transcript *transcript, IDType geneId) { 
  char qStr[1024];
  StatementHandle *sth;
  int i;
  Translation *translation;
  IDType transcriptId;
  IDType xrefId;
  TranslationAdaptor *translAdaptor;
  ExonAdaptor *ea;
  int rank;
  int exonCount;
  

/* NIY
   if( ! ref $transcript || !$transcript->isa('Bio::EnsEMBL::Transcript') ) {
       $self->throw("$transcript is not a EnsEMBL transcript - not dumping!");
   }
*/

   // store translation 
   // then store the transcript 
   // then store the exon_transcript table

  translation = Transcript_getTranslation(transcript);

  ea = DBAdaptor_getExonAdaptor(ta->dba);
  translAdaptor = DBAdaptor_getTranslationAdaptor(ta->dba);

  exonCount = Transcript_getExonCount(transcript);

  for (i=0;i<Transcript_getExonCount(transcript);i++) {
    Exon *exon = Transcript_getExonAt(transcript,i);
    ExonAdaptor_store(ea,exon);
  }

  if (translation) {
    TranslationAdaptor_store(translAdaptor, translation);
  }


  // assuming that the store is used during the Genebuil process, set
  // the display_xref_id to 0.  This ought to get re-set during the protein
  // pipeline run.  This probably update to the gene table has yet to be
  // implemented.
  xrefId = 0;

  // ok - now load this line in
  if (translation) {
    sprintf(qStr,
       "insert into transcript ( gene_id, translation_id, exon_count, display_xref_id )"
       " values ( " IDFMTSTR ", " IDFMTSTR ", %d, 0)",geneId, Translation_getDbID(translation),exonCount);
  } else {
    sprintf(qStr,
       "insert into transcript ( gene_id, translation_id, exon_count, display_xref_id)"
       " values ( " IDFMTSTR ", 0, %d, 0)",geneId, exonCount);
  }
  sth = ta->prepare((BaseAdaptor *)ta,qStr,strlen(qStr));
        
  sth->execute(sth);

  transcriptId = sth->getInsertId(sth);
  sth->finish(sth);

/* NIY
   #print STDERR "Going to look at gene links\n";
   my $dbEntryAdaptor = $self->db->get_DBEntryAdaptor();

   foreach my $dbl ( @{$transcript->get_all_DBLinks} ) {
     $dbEntryAdaptor->store( $dbl, $transc_dbID, "Transcript" );
   }
*/

   // 
   // Update transcript to point to display xref if it is set
   // 
/* NIY
   if(my $dxref = $transcript->display_xref) {
     if(my $dxref_id = $dbEntryAdaptor->exists($dxref)) {
       $self->prepare( "update transcript set display_xref_id = ".
                   $dxref_id . " where transcript_id = ".$transc_dbID);
       $self->execute();
       $dxref->dbID($dxref_id);
       $dxref->adaptor($dbEntryAdaptor);
     }
   }
*/


  sprintf(qStr, "insert into exon_transcript (exon_id,transcript_id,rank) values (%" IDFMTSTR ",%" IDFMTSTR ",%%d)");
  sth = ta->prepare((BaseAdaptor *)ta,qStr,strlen(qStr));
  rank = 1;
  for (i=0; i<exonCount; i++) {
    Exon *exon = Transcript_getExonAt(transcript,i);
    sth->execute(sth, (IDType)Exon_getDbID(exon), (IDType)transcriptId, rank);
    rank++;
  }
  sth->finish(sth);

  if (Transcript_getStableId(transcript)) {
    if (Transcript_getVersion(transcript) == -1) {
      fprintf(stderr, "Error: Trying to store incomplete stable id information for transcript\n");
    }

    sprintf(qStr, "INSERT INTO transcript_stable_id(transcript_id," 
                      "stable_id,version)"
                      " VALUES(" IDFMTSTR ", '%s', %d)",
                      transcriptId,
                      Transcript_getStableId(transcript),
                      Transcript_getVersion(transcript));
    sth = ta->prepare((BaseAdaptor *)ta, qStr, strlen(qStr));
    sth->execute(sth);
    sth->finish(sth);
  }

  Transcript_setDbID(transcript, transcriptId);
  Transcript_setAdaptor(transcript, (BaseAdaptor *)ta);

  return transcriptId;
}

