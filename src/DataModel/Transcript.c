#include "Transcript.h"
#include "TranscriptAdaptor.h"

Transcript *Transcript_new() {
  Transcript *transcript;

  if ((transcript = (Transcript *)calloc(1,sizeof(Transcript))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for transcripe\n");
    return NULL;
  }

  Transcript_setVersion(transcript,-1);

  return transcript;
}

Set *Transcript_getAllDBLinks(Transcript *t) {
  if (!t->dbLinks) {
    TranscriptAdaptor *ta = Transcript_getAdaptor(t);

    if (ta) {
      DBEntryAdaptor *dbea = DBAdaptor_getDBEntryAdaptor(ta->dba);
      t->dbLinks = DBEntryAdaptor_fetchAllByTranscript(dbea,t);
    } else {
      t->dbLinks = emptySet;
    }
  }

  return t->dbLinks;
}

char *Transcript_setType(Transcript *t, char *type) {
  if ((t->type = (char *)malloc(strlen(type)+1)) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for transcript type\n");
    return NULL;
  }

  strcpy(t->type,type);

  return t->type;
}

char *Transcript_getStableId(Transcript *transcript) {
  TranscriptAdaptor *ta = (TranscriptAdaptor *)Transcript_getAdaptor(transcript);

  if (StableIdInfo_getStableId(&(transcript->si)) == NULL && ta) {
    TranscriptAdaptor_getStableEntryInfo(ta,transcript);
  }
  return StableIdInfo_getStableId(&(transcript->si));
}

int Transcript_getVersion(Transcript *transcript) {
  TranscriptAdaptor *ta = (TranscriptAdaptor *)Transcript_getAdaptor(transcript);

  if (StableIdInfo_getVersion(&(transcript->si)) == -1 && ta) {
    TranscriptAdaptor_getStableEntryInfo(ta,transcript);
  }
  return StableIdInfo_getVersion(&(transcript->si));
}

Transcript *Transcript_transform(Transcript *trans, IDHash *exonTransforms) {
  int i;
  Set *mappedExonSet = Set_new();
  
  for (i=0;i<Transcript_getExonCount(trans);i++) {
    Exon *exon = (Exon *)Transcript_getExonAt(trans,i);
    IDType exonRef = (int)exon;

    // the old exon was successfully remapped then store the new exon
/* CHECK */
    if ( IDHash_contains(exonTransforms,exonRef)) {
      Set_addElement(mappedExonSet,IDHash_getValue(exonTransforms, exonRef));
    }
    // but for the case where the exon was unable to be mapped, as it
    // was outside the bounds of the slice, include the original exon.
    else {
      Set_addElement(mappedExonSet,exon);
    }
  }

  //Flush the exons and NIY all related internal caches
  Transcript_flushExons(trans);

  // attach the new list of exons to the transcript
  for (i=0; i<Set_getNumElement(mappedExonSet); i++) {
    Transcript_addExon(trans,(Exon *)Set_getElementAt(mappedExonSet,i));
  }

  if ( Transcript_getTranslation(trans)) {
    Translation_transform(Transcript_getTranslation(trans), exonTransforms);
  }

  Set_free(mappedExonSet,NULL);

  return trans;
}

void Transcript_flushExons(Transcript *trans) {
  Transcript_removeAllExons(trans);
// NIY caches 
}
