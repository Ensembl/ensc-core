#define __TRANSCRIPT_MAIN__
#include "Transcript.h"
#undef __TRANSCRIPT_MAIN__

#include "DBAdaptor.h"
#include "TranscriptAdaptor.h"
#include "TranslationAdaptor.h"
#include "IntronSupportingEvidenceAdaptor.h"
#include "IntronSupportingEvidence.h"
#include "TranscriptSupportingFeatureAdaptor.h"
#include "AttributeAdaptor.h"
#include "DBEntryAdaptor.h"
#include "StrUtil.h"
#include "SeqUtil.h"
#include "translate.h"

#include "Attribute.h"

#include <stdlib.h>

Transcript *Transcript_new() {
  Transcript *transcript;

  if ((transcript = (Transcript *)calloc(1,sizeof(Transcript))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for transcripe\n");
    return NULL;
  }

  Transcript_setModified(transcript,0);
  Transcript_setCreated(transcript,0);
  Transcript_setVersion(transcript,-1);
  Transcript_setIsCurrent(transcript,1);
// -1 is invalid is canonical value, when call getIsCanonical value will be retrieved from db
  Transcript_setIsCanonical(transcript,-1);

  transcript->exons = Vector_new();

  transcript->objectType = CLASS_TRANSCRIPT;
  Object_incRefCount(transcript);

  transcript->funcs = &transcriptFuncs;

  return transcript;
}

Transcript *Transcript_shallowCopy(Transcript *transcript) {
  Transcript *newTranscript = Transcript_new();

  memcpy(newTranscript,transcript,sizeof(Transcript));

  return newTranscript;
}

/*
=head2 swap_exons

  Arg [1]    : Bio::EnsEMBL::Exon $old_Exon
               An exon that should be replaced
  Arg [2]    : Bio::EnsEMBL::Exon $new_Exon
               The replacement Exon
  Example    : none
  Description: exchange an exon in the current Exon list with a given one.
               Usually done before storing of Gene, so the Exons can
               be shared between Transcripts.
  Returntype : none
  Exceptions : none
  Caller     : GeneAdaptor->store()
  Status     : Stable

=cut
*/
// New
// I tried to match the internal v external access to variables
void Transcript_swapExons(Transcript *transcript, Exon *oldExon, Exon *newExon) {

  int i;
  for (i=0; i<Vector_getNumElement(transcript->exons); i++) {
    Exon *transExon = Vector_getElementAt(transcript->exons, i);
    if (transExon == oldExon) {
      Vector_setElementAt(transcript->exons, i, newExon);
      break;
    }
  }

  if (transcript->translation) {
    if (Translation_getStartExon(transcript->translation) == oldExon) {
      Translation_setStartExon(transcript->translation, newExon);
    }
    if (Translation_getEndExon(transcript->translation) == oldExon) {
      Translation_setEndExon(transcript->translation, newExon);
    }
  }
}


/*
=head2 is_canonical

  Args [1]      : (optional) Boolean is_canonical

  Example       : if ($transcript->is_canonical()) { ... }

  Description : Returns true (non-zero) if the transcript is the
                canonical transcript of its gene, false (0) if not. If the code
                returns an undefined it is because its state is not currently
                known. Internally the code will consult the database for this
                value if it is unknown and the transcript has a dbID and an
                attached adaptor

  Return type   : Boolean

  Status        : Stable

=cut
*/
// New
int Transcript_getIsCanonical(Transcript *transcript) {
  //Shortcut call
  if (transcript->isCanonical >= 0) {
    return transcript->isCanonical;
  }

  if (Transcript_getDbID(transcript) && Transcript_getAdaptor(transcript) != NULL) {
    
    TranscriptAdaptor *ta = (TranscriptAdaptor *)Transcript_getAdaptor(transcript);
    transcript->isCanonical = TranscriptAdaptor_isTranscriptCanonical(ta, transcript);
  } else {
    fprintf(stderr, "Called getIsCanonical when isCanonical is undefined and don't have an adaptor. Can't determine whether this is canonical - bye!\n");
    exit(1);
  }

  return transcript->isCanonical;
}

/*
=head2 get_all_Attributes

  Arg [1]    : optional string $attrib_code
               The code of the attribute type to retrieve values for.
  Example    : ($rna_edits) = @{$transcript->get_all_Attributes('_rna_edit')};
               @transc_attributes    = @{$transcript->get_all_Attributes()};
  Description: Gets a list of Attributes of this transcript.
               Optionally just get Attrubutes for given code.
  Returntype : listref Bio::EnsEMBL::Attribute
  Exceptions : warning if transcript does not have attached adaptor and
               attempts lazy load.
  Caller     : general
  Status     : Stable

=cut
*/
// New

// NIY:
// Because this can filter the results the vector that gets returned must be freeable - so for now
// make a copy of the translation->attributes vector if returning unfiltered so behaviour is
// consistent. Long term probably want reference count incremented
Vector *Transcript_getAllAttributes(Transcript *transcript, char *attribCode) {
  if (transcript->attributes == NULL) {
    TranscriptAdaptor *ta = (TranscriptAdaptor *)Transcript_getAdaptor(transcript);
    if (ta == NULL) { // No adaptor
// Perl comments out the warning, I'll put it back for now, just in case
      fprintf(stderr,"Warning: Cannot get attributes without an adaptor.\n");
      return Vector_new();
    }

    AttributeAdaptor *ata = DBAdaptor_getAttributeAdaptor(ta->dba);
    transcript->attributes = AttributeAdaptor_fetchAllByTranscript(ata, transcript, NULL);
  }

  if (attribCode != NULL) {
    Vector *results = Vector_new();
    int i;
    for (i=0; i<Vector_getNumElement(transcript->attributes); i++) {
      Attribute *attrib = Vector_getElementAt(transcript->attributes, i);
      if (!strcasecmp(attrib->code, attribCode)) {
        Vector_addElement(results, attrib);
      }
    }
    return results;
  } else {
// See NIY note above for why I'm making a copy
    return Vector_copy(transcript->attributes);
  }
}


// Note this used to be called getAllDBLinks but it looks more like the current getAllDBEntries so call it that
Vector *Transcript_getAllDBEntries(Transcript *t) {
  if (!t->dbLinks) {
    TranscriptAdaptor *ta = (TranscriptAdaptor *)Transcript_getAdaptor(t);

    if (ta) {
      DBEntryAdaptor *dbea = DBAdaptor_getDBEntryAdaptor(ta->dba);
      DBEntryAdaptor_fetchAllByTranscript(dbea,t);
    } else {
      t->dbLinks = emptyVector;
    }
  }

  return t->dbLinks;
}

int Transcript_addDBLink(Transcript *t, DBEntry *dbe) {
  if (!t->dbLinks) {
    t->dbLinks = Vector_new();
  }

  Vector_addElement(t->dbLinks, dbe); 
  return 1;
}

ECOSTRING Transcript_setBiotype(Transcript *t, char *biotype) {
  EcoString_copyStr(ecoSTable, &(t->biotype),biotype,0);

  if (t->biotype == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for biotype\n");
    return NULL;
  }

  return t->biotype;
}

ECOSTRING Transcript_setStatus(Transcript *t, char *status) {
  EcoString_copyStr(ecoSTable, &(t->status),status,0);

  if (t->status == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for status\n");
    return NULL;
  }

  return t->status;
}

char *Transcript_setDescription(Transcript *t, char *description) {
  if ((t->description = (char *)malloc(strlen(description)+1)) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for description\n");
    return NULL;
  }

  strcpy(t->description,description);

  return t->description;
}

ECOSTRING Transcript_setExternalDb(Transcript *t, char *externalDb) {
  EcoString_copyStr(ecoSTable, &(t->externalDb),externalDb,0);

  if (t->externalDb == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for externalDb\n");
    return NULL;
  }

  return t->externalDb;
}

ECOSTRING Transcript_setExternalStatus(Transcript *t, char *externalStatus) {
  EcoString_copyStr(ecoSTable, &(t->externalStatus),externalStatus,0);

  if (t->externalStatus == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for externalStatus\n");
    return NULL;
  }

  return t->externalStatus;
}

char *Transcript_setExternalName(Transcript *t, char *externalName) {
  if (externalName == NULL) {
    t->externalName = NULL;
    return NULL;
  }
  if ((t->externalName = (char *)malloc(strlen(externalName)+1)) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for externalName\n");
    return NULL;
  }

  strcpy(t->externalName,externalName);

  return t->externalName;
}

char *Transcript_getStableId(Transcript *transcript) {
  TranscriptAdaptor *ta = (TranscriptAdaptor *)Transcript_getAdaptor(transcript);

  if (StableIdInfo_getStableId(&(transcript->si)) == NULL && ta) {
//    TranscriptAdaptor_getStableEntryInfo(ta,transcript);
  }
  return StableIdInfo_getStableId(&(transcript->si));
}

int Transcript_getVersion(Transcript *transcript) {
  TranscriptAdaptor *ta = (TranscriptAdaptor *)Transcript_getAdaptor(transcript);

  if (StableIdInfo_getVersion(&(transcript->si)) == -1 && ta) {
//    TranscriptAdaptor_getStableEntryInfo(ta,transcript);
  }
  return StableIdInfo_getVersion(&(transcript->si));
}

/*
=head2 transfer

  Arg  1     : Bio::EnsEMBL::Slice $destination_slice
  Example    : $transcript = $transcript->transfer($slice);
  Description: Moves this transcript to the given slice.
               If this Transcripts has Exons attached, they move as well.
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/
// New
int notWarned = 1;
Transcript *Transcript_transfer(Transcript *transcript, Slice *slice) {
  // Call super transfer
  Transcript *newTranscript = SeqFeature_transfer(transcript, slice);

  if (newTranscript == NULL) {
    return NULL;
  }

  if (transcript->translation) {
    Translation *newTranslation = Translation_new();
    memcpy(newTranslation, transcript->translation, sizeof(Translation));
//
//    Transcript_setTranslation(new_transcript->{'translation'} = $new_translation;
    // Perl did direct set so maybe I'll do same
    newTranscript->translation = newTranslation;
  }

//
//Perl  if ( exists $self->{'_trans_exon_array'} ) {
  if (transcript->exons != NULL && Vector_getNumElement(transcript->exons)) {
    Vector *newExons = Vector_new();

    int i;
    for (i=0;i<Vector_getNumElement(transcript->exons);i++) {
      Exon *oldExon = Vector_getElementAt(transcript->exons, i);
      Exon *newExon = Exon_transfer(oldExon, slice);

      if (newTranscript->translation) {
        Translation *newTranslation = Transcript_getTranslation(newTranscript);

        if( Translation_getStartExon(newTranslation) == oldExon ) {
          Translation_setStartExon(newTranslation, newExon);
        }
        if( Translation_getEndExon(newTranslation) == oldExon ) {
          Translation_setEndExon(newTranslation, newExon);
        }
      }
      Vector_addElement(newExons, newExon );
    }

    newTranscript->exons = newExons;

// NIY Free old stuff

  }

  if (notWarned) {
    fprintf(stderr,"cache clearing not implemented yet for transcript transfer\n");
    notWarned = 0;
  }

  if (transcript->supportingEvidence != NULL && Vector_getNumElement(transcript->supportingEvidence) != 0) {
    Vector *newFeatures = Vector_new();

    int i;
    for (i=0; i<Vector_getNumElement(transcript->supportingEvidence); i++) {
      SeqFeature *oldFeature = Vector_getElementAt(transcript->supportingEvidence, i);
      SeqFeature *newFeature = SeqFeature_transfer(oldFeature, slice);
      Vector_addElement(newFeatures, newFeature);
    }
    newTranscript->supportingEvidence = newFeatures;
  }

  if (transcript->iseVector != NULL && Vector_getNumElement(transcript->iseVector)) {
    Vector *newFeatures = Vector_new();

    int i;
    for (i=0; i<Vector_getNumElement(transcript->iseVector); i++) {
      SeqFeature *oldFeature = Vector_getElementAt(transcript->iseVector, i);
      SeqFeature *newFeature = SeqFeature_transfer(oldFeature, slice);
      Vector_addElement(newFeatures, newFeature);
    }
    newTranscript->iseVector = newFeatures;
  }


/*
  // flush cached internal values that depend on the exon coords
  $new_transcript->{'transcript_mapper'}   = undef;
  $new_transcript->{'coding_region_start'} = undef;
  $new_transcript->{'coding_region_end'}   = undef;
  $new_transcript->{'cdna_coding_start'}   = undef;
  $new_transcript->{'cdna_coding_end'}     = undef;
*/

  return newTranscript;
}

//NIY Freeing old exons???
Transcript *Transcript_transform(Transcript *trans, IDHash *exonTransforms) {
  Vector *mappedExonVector = Vector_new();
  int i;
  
  for (i=0;i<Transcript_getExonCount(trans);i++) {
    Exon *exon = (Exon *)Transcript_getExonAt(trans,i);
    IDType exonRef = (IDType)exon;

    // the old exon was successfully remapped then store the new exon
/* CHECK */
    if ( IDHash_contains(exonTransforms,exonRef)) {
      Vector_addElement(mappedExonVector,IDHash_getValue(exonTransforms, exonRef));
    }
    // but for the case where the exon was unable to be mapped, as it
    // was outside the bounds of the slice, include the original exon.
    else {
      Vector_addElement(mappedExonVector,exon);
    }
  }

  //Flush the exons and NIY all related internal caches
  Transcript_flushExons(trans);

  // attach the new list of exons to the transcript
  for (i=0; i<Vector_getNumElement(mappedExonVector); i++) {
// Added rank arg temporarily just to get to compile
    Transcript_addExon(trans,(Exon *)Vector_getElementAt(mappedExonVector,i),0);
  }

  if ( Transcript_getTranslation(trans)) {
    Translation_transform(Transcript_getTranslation(trans), exonTransforms);
  }

  Vector_free(mappedExonVector);

  return trans;
}

/*
=head2 get_all_supporting_features

  Example    : my @evidence = @{ $transcript->get_all_supporting_features };
  Description: Retreives any supporting features added manually by
               calls to add_supporting_features.
  Returntype : Listref of Bio::EnsEMBL::FeaturePair objects
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/
// New
Vector *Transcript_getAllSupportingFeatures(Transcript *transcript) {

  if( transcript->supportingEvidence == NULL && Transcript_getAdaptor(transcript) != NULL) {
    TranscriptAdaptor *adaptor = (TranscriptAdaptor *)Transcript_getAdaptor(transcript);
    TranscriptSupportingFeatureAdaptor *tsfa = DBAdaptor_getTranscriptSupportingFeatureAdaptor(adaptor->dba);
 
    transcript->supportingEvidence = TranscriptSupportingFeatureAdaptor_fetchAllByTranscript(tsfa, transcript);
  }
  //return $self->{_supporting_evidence} || [];
  return transcript->supportingEvidence;
}



/*
=head2 get_all_IntronSupportingEvidence

  Example     : $ise->get_all_IntronSupportingEvidence();
  Description : Fetches all ISE instances linked to this Transript
  Returntype  : ArrayRef[Bio::EnsEMBL::IntronSupportEvidence] retrieved from
                the DB or from those added via C<add_IntronSupportingEvidence>
  Exceptions  : None

=cut
*/
// New
Vector *Transcript_getAllIntronSupportingEvidence(Transcript *transcript) {

  if( transcript->iseVector == NULL && Transcript_getAdaptor(transcript) != NULL) {
    TranscriptAdaptor *adaptor = (TranscriptAdaptor *)Transcript_getAdaptor(transcript);
    IntronSupportingEvidenceAdaptor *isea = DBAdaptor_getIntronSupportingEvidenceAdaptor(adaptor->dba);
 
    transcript->iseVector = IntronSupportingEvidenceAdaptor_fetchAllByTranscript(isea, transcript);
  }
  return transcript->iseVector;
}

/*
=head2 add_IntronSupportingEvidence

  Arg [1]     : Bio::EnsEMBL::IntronSupportEvidence Object to add
  Example     : $ise->add_IntronSupportingEvidence($ise);
  Description : Adds the IntronSupportEvidence instance to this Transcript. The
                code checks to see if it is a unique ISE instance
  Returntype  : Boolean; true means it was added. False means it was not
                as this ISE was already attached
  Exceptions  : None

=cut
*/
// New, but half hearted equals check implementation
int Transcript_addIntronSupportingEvidence(Transcript *transcript, IntronSupportingEvidence *ise) {
  int unique = 1;

  if (transcript->iseVector == NULL) {
    transcript->iseVector = Vector_new();
  } else {
    int i;
    for (i=0; i<Vector_getNumElement(transcript->iseVector); i++) {
      IntronSupportingEvidence *compIse = Vector_getElementAt(transcript->iseVector, i);
      
//Equals check not same as perl - I don't care that much and its acres of code to do these checks
      if ( compIse == ise || 
           ( IntronSupportingEvidence_getStart(ise) == IntronSupportingEvidence_getStart(ise) &&
             IntronSupportingEvidence_getEnd(ise) == IntronSupportingEvidence_getEnd(ise) &&
             IntronSupportingEvidence_getStrand(ise) == IntronSupportingEvidence_getStrand(ise) &&
             Slice_getSeqRegionId(IntronSupportingEvidence_getSlice(ise)) == Slice_getSeqRegionId(IntronSupportingEvidence_getSlice(ise)) &&
             IntronSupportingEvidence_getAnalysis(ise) == IntronSupportingEvidence_getAnalysis(ise)
           )
         ) {
        unique = 0;
        break;
      }
    }
  }

  if (unique) {
    Vector_addElement(transcript->iseVector, ise);
    return 1;
  }
  return 0;
}


void Transcript_flushExons(Transcript *trans) {
  Transcript_removeAllExons(trans);
// NIY caches 
}


int Transcript_setCodingRegionEnd(Transcript *trans, int end) {
  trans->codingRegionEnd = end;
  Transcript_setCodingRegionEndIsSet(trans,TRUE);
  return trans->codingRegionEnd;
}

int Transcript_getCodingRegionEnd(Transcript *trans) {
  Translation *translation = Transcript_getTranslation(trans);

  if (!translation) {
    return POS_UNDEF;
  }

  if (!Transcript_getCodingRegionEndIsSet(trans) && translation) {
    int end;
    int strand = Exon_getStrand(Translation_getStartExon(translation));

    if (strand == 1) {
      end = Exon_getStart(Translation_getEndExon(translation));
      end += (Translation_getEnd(translation) - 1);
    } else {
      end = Exon_getEnd(Translation_getStartExon(translation));
      end -= (Translation_getStart(translation) - 1 );
    }
    Transcript_setCodingRegionEnd(trans,end);
  }

  return trans->codingRegionEnd;
}

int Transcript_setCodingRegionStart(Transcript *trans, int start) {
  trans->codingRegionStart = start;
  Transcript_setCodingRegionStartIsSet(trans,TRUE);
  return trans->codingRegionStart;
}

int Transcript_getCodingRegionStart(Transcript *trans) {
  Translation *translation = Transcript_getTranslation(trans);

  if (!translation) {
    return POS_UNDEF;
  }

  if (!Transcript_getCodingRegionStartIsSet(trans) && translation) {
    int start;
    int strand = Exon_getStrand(Translation_getStartExon(translation));

    if (strand == 1) {
      start = Exon_getStart(Translation_getStartExon(translation));
      start += (Translation_getStart(translation) - 1);
    } else {
      start = Exon_getEnd(Translation_getEndExon(translation));
      start -= (Translation_getEnd(translation) - 1 );
    }
    Transcript_setCodingRegionStart(trans,start);
  }

  return trans->codingRegionStart;
}

void Transcript_sort(Transcript *trans) {
  Exon *firstExon = Transcript_getExonAt(trans,0);
  int strand = Exon_getStrand(firstExon);

  if (strand == 1) {
    Transcript_sortExons(trans, Exon_forwardStrandCompFunc);
  } else if (strand == -1) {
    Transcript_sortExons(trans, Exon_reverseStrandCompFunc);
  }
}

Exon *Transcript_getStartExon(Transcript *trans) {
  return Transcript_getExonAt(trans,0); 
}

Exon *Transcript_getEndExon(Transcript *trans) {
  return Transcript_getExonAt(trans,Transcript_getExonCount(trans)-1); 
}

Translation *Transcript_getTranslation(Transcript *trans) {
  if (!trans->translation && trans->translationId ) {
    TranslationAdaptor *ta = DBAdaptor_getTranslationAdaptor(Transcript_getAdaptor(trans)->dba);
    
    trans->translation = TranslationAdaptor_fetchByTranscript(ta, trans);
  }
  return trans->translation;
}


char *Transcript_getSplicedSeq(Transcript *trans) {
  int i;
  char *seqStr;

  
  StrUtil_copyString(&seqStr, "", 0);

  for (i=0; i<Transcript_getExonCount(trans); i++) {
    Exon *ex = Transcript_getExonAt(trans, i);
    seqStr = StrUtil_appendString(seqStr,Exon_getSeqString(ex));
  }


  return seqStr;
}

char *Transcript_getTranslateableSeq(Transcript *trans) {
  int i;
  int lastPhase = 0;
  int first = 1;
  char *mrna = StrUtil_copyString(&mrna,"",0);
  Vector *translateableExons;

  translateableExons = Transcript_getAllTranslateableExons(trans);
 
  for (i=0; i<Vector_getNumElement(translateableExons); i++) {
    Exon *exon = Vector_getElementAt(translateableExons, i);
    //fprintf(stderr, "translateable exons %d from %ld to %ld\n", i, Exon_getStart(exon), Exon_getEnd(exon));
    int phase = 0;

// NIY    if (defined($exon->phase)) {
      phase = Exon_getPhase(exon);
//    }

    if (first) {
      mrna = SeqUtil_addNs(mrna,phase);
      first = 0;
    }

/*
    if (phase != lastPhase && getenv("MONKEY_EXONS")) {
      // endpadding for the last exon
      if (lastPhase == 1 ) {
        mrna = StrUtil_appendString(mrna,"NN");
      } else if (lastPhase == 2) {
        mrna = StrUtil_appendString(mrna,"N");
      }
      //startpadding for this exon
      mrna = SeqUtil_addNs(mrna,phase);
    }
*/
    mrna = StrUtil_appendString(mrna, Exon_getSeqString(exon));
    lastPhase = Exon_getEndPhase(exon);
// Hack to test memory usage
//    free(Exon_getSeqCacheString(exon));
//    Exon_setSeqCacheString(exon, NULL);
  }

  Vector_free(translateableExons);

  //fprintf(stderr,"Length of seq = %d\n", strlen(mrna));

// NIY free adjusted

  return mrna;
}


/*
=head2 cdna_coding_start

  Arg [1]    : (optional) $value
  Example    : $relative_coding_start = $transcript->cdna_coding_start;
  Description: Retrieves the position of the coding start of this transcript
               in cdna coordinates (relative to the start of the 5prime end of
               the transcript, excluding introns, including utrs).
  Returntype : int
  Exceptions : none
  Caller     : five_prime_utr, get_all_snps, general

=cut
*/

int Transcript_setcDNACodingStart(Transcript *trans, int start) {
  trans->cDNACodingStart = start;
  Transcript_setcDNACodingStartIsSet(trans,TRUE);
  return trans->cDNACodingStart;
}

int Transcript_getcDNACodingStart(Transcript *trans) {
  Translation *translation = Transcript_getTranslation(trans);

  if (!Transcript_getcDNACodingStartIsSet(trans) && translation) {
    int start = 0;
    int i;

    for (i=0; i<Transcript_getExonCount(trans); i++) {
      Exon *exon = Transcript_getExonAt(trans,i);

      if (exon == Translation_getStartExon(translation)) {
        //add the utr portion of the start exon
        start += Translation_getStart(translation);
        break;
      } else {
        //add the entire length of this non-coding exon
        start += Exon_getLength(exon);
      }
    }
    Transcript_setcDNACodingStart(trans,start);
  }

  return trans->cDNACodingStart;
}

/*
=head2 cdna_coding_end

  Arg [1]    : (optional) $value
  Example    : $cdna_coding_end = $transcript->cdna_coding_end;
  Description: Retrieves the end of the coding region of this transcript in
               cdna coordinates (relative to the five prime end of the
               transcript, excluding introns, including utrs).
               Note
  Returntype : none
  Exceptions : none
  Caller     : general

=cut
*/

int Transcript_setcDNACodingEnd(Transcript *trans, int end) {
  trans->cDNACodingEnd = end;
  Transcript_setcDNACodingEndIsSet(trans,TRUE);
  return trans->cDNACodingEnd;
}

int Transcript_getcDNACodingEnd(Transcript *trans) {
  Translation *translation = Transcript_getTranslation(trans);

  if (!Transcript_getcDNACodingEndIsSet(trans) && translation) {
    int end = 0;
    int i;

    for (i=0; i<Transcript_getExonCount(trans); i++) {
      Exon *exon = Transcript_getExonAt(trans,i);

      if (exon == Translation_getEndExon(translation)) {
        //add the coding portion of the final coding exon
        end += Translation_getEnd(translation);
        break;
      } else {
        //add the entire length
        end += Exon_getLength(exon);
      }
    }
    Transcript_setcDNACodingEnd(trans,end);
  }

  return trans->cDNACodingEnd;
}

int Transcript_getLength(Transcript *trans) {
  int length = 0;
  int i;

  for (i=0; i<Transcript_getExonCount(trans); i++) {
    Exon *ex = Transcript_getExonAt(trans, i);
    length += Exon_getLength(ex);
  }
  return length;
}


Vector *Transcript_getAllTranslateableExons(Transcript *trans) {
  Vector *translateable;
  Translation *translation;
  Exon *startExon;
  Exon *endExon;
  int tlnStart;
  int tlnEnd;
  int i;
  

  if (!(translation = Transcript_getTranslation(trans))) {
    fprintf(stderr, "Error: No translation attached to transcript object\n");
    exit(1);
  }
  startExon      = Translation_getStartExon(translation);
  endExon        = Translation_getEndExon(translation);
  tlnStart       = Translation_getStart(translation);
  tlnEnd         = Translation_getEnd(translation);

// Hack
  Transcript_sort(trans);

  translateable = Vector_new();

  //fprintf(stderr,"Looking for start exon %p %s\n",startExon, Exon_getStableId(startExon));
  //fprintf(stderr,"Exon count %d\n",Transcript_getExonCount(trans));
  for (i=0; i<Transcript_getExonCount(trans); i++) {
    Exon *ex = Transcript_getExonAt(trans,i);
    int length;
    int adjustStart = 0;
    int adjustEnd   = 0;

    //fprintf(stderr,"   exon %p %s\n",ex, Exon_getStableId(ex));
 
    if (ex != startExon && !Vector_getNumElement(translateable)) {
       //fprintf(stderr,"  skip %p %s\n",ex, Exon_getStableId(ex));
 
      continue;   // Not yet in translated region
    }

    length  = Exon_getLength(ex);

    // Adjust to translation start if this is the start exon
    if (ex == startExon ) {
      if (tlnStart < 1 || tlnStart > length) {
        fprintf(stderr,"Error: Translation start %d is outside exon %p length=%d\n", tlnStart,ex,length);
        exit(1);
      }
      adjustStart = tlnStart - 1;
    }

    // Adjust to translation end if this is the end exon
    if (ex == endExon) {
      if (tlnEnd < 1 || tlnEnd > length) {
        fprintf(stderr,"Error: Translation start %d is outside exon %p length=%d\n", tlnEnd,ex,length);
        exit(1);
      }
      adjustEnd = tlnEnd - length;
    }

    // Make a truncated exon if the translation start or
    // end causes the coordinates to be altered.
    if (adjustEnd || adjustStart) {
      Exon *newex = Exon_adjustStartEnd(ex, adjustStart, adjustEnd);

      Vector_addElement(translateable, newex );
    } else {
      Vector_addElement(translateable, ex);
    }

    // Exit the loop when we've found the last exon
    if (ex == endExon) break;
  }
  //fprintf(stderr," Have %d translateable exons\n", Vector_getNumElement(translateable));
  return translateable;
}


char *Transcript_translate(Transcript *trans) {
  char *mRNA;
  int lenmRNA;
  int lengths[6];
  char *frm[6];
  int i;

  mRNA = Transcript_getTranslateableSeq(trans);

  long codonTableId = 1; // Default to 'standard' table

  Slice *slice = Transcript_getSlice(trans);
  if (slice) {
    Vector *ctAttribs = Slice_getAllAttributes(slice, "codon_table");
    // NIY free properly - Slice attributes current are NOT cached in Slice
    if (Vector_getNumElement(ctAttribs) > 1) {
      fprintf(stderr,"Invalid number of codon_table attributes for slice %s\n", Slice_getName(slice));
      exit(1);
    }
    if (Vector_getNumElement(ctAttribs) == 1) {
      Attribute *attrib = Vector_getElementAt(ctAttribs, 0);
      if (!StrUtil_isLongInteger(&codonTableId, Attribute_getValue(attrib))) {
        fprintf(stderr,"Invalid codon_table attribute format for slice %s attribute value %s\n", Slice_getName(slice),Attribute_getValue(attrib));
        exit(1);
      }

      Vector_setFreeFunc(ctAttribs, Attribute_free);
    }
    Vector_free(ctAttribs);
  }
  
    
  lenmRNA = strlen(mRNA);

  if (lenmRNA%3 == 0 && lenmRNA) {
    char *lastCodon = &(mRNA[lenmRNA-3]);
    if (!strcmp(lastCodon,"TAG") ||
        !strcmp(lastCodon,"TGA") ||
        !strcmp(lastCodon,"TAA")) {
      *lastCodon = '\0';
    }
  }
  // the above lines will remove the final stop codon from the mRNA
  // sequence produced if it is present, this is so any peptide produced
  // won't have a terminal stop codon
  // if you want to have a terminal stop codon either comment this line out
  // or call translatable seq directly and produce a translation from it

  for (i=0;i<6;i++) {
    frm[i] = (char *)malloc(lenmRNA/3 + 2);
  }

  //fprintf(stderr, "translateable seq = %s\n",mRNA);
  translate(mRNA, frm, lengths, codonTableId);

  free(mRNA);

  // Note don't free frame 0
  for (i=1;i<6;i++) {
    free(frm[i]);
  }

  return Translation_modifyTranslation(Transcript_getTranslation(trans), frm[0]);
}


#ifdef DONE
char *Transcript_seq(Transcript *trans) {
  int i;
  char *transcriptSeqString;

  for (i=0; i<Transcript_getExonCount(trans); i++) {
    Exon *ex = Transcript_getExonAt(trans,i);
    transcriptSeqString = StrUtil_appendString(transcriptSeqString, Exon_getSeqString(ex));
  }

    my $seq = Bio::Seq->new(
        -DISPLAY_ID => $self->stable_id,
        -MOLTYPE    => 'dna',
        -SEQ        => $transcript_seq_string,
        );

  return seq;
}
#endif

MapperRangeSet *Transcript_pep2Genomic(Transcript *trans, int start, int end) {

  // move start end into translate cDNA coordinates now.
  // much easier!
  start = 3*start-2 + (Transcript_getcDNACodingStart(trans) - 1);
  end   = 3*end + (Transcript_getcDNACodingStart(trans) - 1);

  return Transcript_cDNA2Genomic(trans, start, end);
}


MapperRangeSet *Transcript_genomic2Pep(Transcript *trans, int start, int end, int strand, BaseContig *contig) {
  MapperRangeSet *mapped;
  MapperRangeSet *out;
  int startPhase;
  int i;

  mapped = Transcript_genomic2cDNA(trans, start, end, strand, contig);

  out = MapperRangeSet_new();

  if(Transcript_getExonCount(trans)) {
    Exon *exon = Transcript_getExonAt(trans,0);
    startPhase = Exon_getPhase(exon);
  } else {
    startPhase = -1;
  }

  for (i=0;i<mapped->nRange;i++) {
    MapperRange *mr = MapperRangeSet_getRangeAt(mapped,i);
    if (mr->rangeType == MAPPERRANGE_GAP) {
      MapperRangeSet_addRange(out, mr);
    } else {
      MapperCoordinate *coord = (MapperCoordinate *)mr;
      int start = coord->start;
      int end   = coord->end;
      int cdnaCStart = Transcript_getcDNACodingStart(trans);
      int cdnaCEnd   = Transcript_getcDNACodingEnd(trans);
      
      if (coord->strand == -1 || end < cdnaCStart || start > cdnaCEnd) {
	// is all gap - does not map to peptide
	MapperGap *gap = MapperGap_new(start,end,0);
        MapperRangeSet_addRange(out, (MapperRange *)gap);
      } else {
	// we know area is at least partially overlapping CDS
	int cdsStart = start - cdnaCStart + 1;
	int cdsEnd   = end   - cdnaCStart + 1;
	MapperGap *endGap = NULL;
        int shift;
        int pepStart;
        int pepEnd;
        MapperCoordinate *newCoord;

	if (start < cdnaCStart) {
	  // start of coordinates are in the 5prime UTR
	  MapperGap *gap = MapperGap_new(start,cdnaCStart-1,0);
	  // start is now relative to start of CDS
	  cdsStart = 1;
          MapperRangeSet_addRange(out, (MapperRange *)gap);
	} 
	
	if (end > cdnaCEnd) {
	  // end of coordinates are in the 3prime UTR
	  endGap = MapperGap_new(cdnaCEnd+1, end,0);
	  // adjust end to relative to CDS start
	  cdsEnd = cdnaCEnd - cdnaCStart + 1;
	}

	// start and end are now entirely in CDS and relative to CDS start

	// take into account possible N padding at beginning of CDS
	shift = (startPhase > 0) ? startPhase : 0;
	
	// convert to peptide coordinates
        // NOTE Don't use same coord unlike perl so can easily free mapped
	pepStart = (cdsStart + shift + 2) / 3;
	pepEnd   = (cdsEnd   + shift + 2) / 3;
  //NIY: What should coordsystem be???
        newCoord = MapperCoordinate_new(coord->id, pepStart, pepEnd, coord->strand,NULL, 0);
	
        MapperRangeSet_addRange(out, (MapperRange *)newCoord);

	if (endGap) {
	  // push out the region which was in the 3prime utr
          MapperRangeSet_addRange(out, (MapperRange *)endGap);
	}
      }	
    }
  }

// NIY Free mapped, but remember bits have been used 
  MapperRangeSet_free(mapped);

  return out;
}

MapperRangeSet *Transcript_cDNA2Genomic(Transcript *trans, int start, int end) {
  Mapper *mapper;

  mapper = Transcript_getcDNACoordMapper(trans);

  return Mapper_mapCoordinates(mapper, (IDType)trans, start, end, 1, "cdna" );
}

MapperRangeSet *Transcript_genomic2cDNA(Transcript *trans, int start, int end, int strand, BaseContig *contig) {
  Mapper *mapper;

  // "ids" in mapper are contigs of exons, so use the same contig that should
  // be attached to all of the exons...
  Exon *firstExon = Transcript_getExonAt(trans,0);
  if (!contig) {
    contig = Exon_getContig(firstExon);
  }

  mapper = Transcript_getcDNACoordMapper(trans);


  printf("MAPPING %d - %d (%d)\n",start,end,strand);
  printf(" %s = %s\n", BaseContig_getName(contig), 
         BaseContig_getName(Exon_getContig(firstExon)));

  return Mapper_mapCoordinates(mapper,(IDType)contig, start, end, strand, "genomic");
}


Mapper *Transcript_getcDNACoordMapper(Transcript *trans) {
  Mapper *mapper;
  int start = 1;
  int i;

  if (trans->exonCoordMapper) {
    return trans->exonCoordMapper;
  }

  // 
  // the mapper is loaded with OBJECTS in place of the IDs !!!!
  //  the objects are the contigs in the exons
  // 
// NIY: What should coordsystems be
  mapper = Mapper_new( "cdna", "genomic", NULL, NULL );

  for (i=0; i<Transcript_getExonCount(trans); i++) {
    Exon *exon = Transcript_getExonAt(trans,i);

    Exon_loadGenomicMapper(exon, mapper, (IDType)trans, start);
    start += Exon_getLength(exon);
  }
  trans->exonCoordMapper = mapper;
  return mapper;
}

void Transcript_free(Transcript *trans) {
  Object_decRefCount(trans);

  if (Object_getRefCount(trans) > 0) {
    return;
  } else if (Object_getRefCount(trans) < 0) {
    fprintf(stderr,"Error: Negative reference count for Transcript\n"
                   "       Freeing it anyway\n");
  }

//  printf("Transcript_free not implemented\n");
}
