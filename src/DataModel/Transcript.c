#include "Transcript.h"
#include "TranscriptAdaptor.h"
#include "DBEntryAdaptor.h"

Transcript *Transcript_new() {
  Transcript *transcript;

  if ((transcript = (Transcript *)calloc(1,sizeof(Transcript))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for transcripe\n");
    return NULL;
  }

  Transcript_setVersion(transcript,-1);

  transcript->objectType = CLASS_TRANSCRIPT;

  return transcript;
}

Vector *Transcript_getAllDBLinks(Transcript *t) {
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
  Vector *mappedExonVector = Vector_new();
  
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
    Transcript_addExon(trans,(Exon *)Vector_getElementAt(mappedExonVector,i));
  }

  if ( Transcript_getTranslation(trans)) {
    Translation_transform(Transcript_getTranslation(trans), exonTransforms);
  }

  Vector_free(mappedExonVector,NULL);

  return trans;
}

void Transcript_flushExons(Transcript *trans) {
  Transcript_removeAllExons(trans);
// NIY caches 
}

int Transcript_setStart(SeqFeature *sf, int start) {
  Transcript *trans = (Transcript *)sf;
  trans->start = start;
  Transcript_setStartIsSet(trans,TRUE);
  return trans->start;
}

int Transcript_getStart(SeqFeature *sf) {
  Transcript *trans = (Transcript *)sf;
  if (!Transcript_getStartIsSet(trans)) {
    int start;
    int strand = Exon_getStrand(Transcript_getStartExon(trans));

    if (strand == 1) {
      start = Exon_getStart(Transcript_getStartExon(trans));
    } else {
      start = Exon_getStart(Transcript_getEndExon(trans));
    }
    trans->start = start;
    Transcript_setStartIsSet(trans,TRUE);
  }

  return trans->start;
}


int Transcript_setEnd(SeqFeature *sf, int end) {
  Transcript *trans = (Transcript *)sf;
  trans->end = end;
  Transcript_setEndIsSet(trans,TRUE);
  return trans->end;
}

int Transcript_getEnd(SeqFeature *sf) {
  Transcript *trans = (Transcript *)sf;
  if (!Transcript_getEndIsSet(trans)) {
    int end;
    int strand = Exon_getStrand(Transcript_getStartExon(trans));

    if (strand == 1) {
      end = Exon_getEnd(Transcript_getEndExon(trans));
    } else {
      end = Exon_getEnd(Transcript_getStartExon(trans));
    }
    trans->end = end;
    Transcript_setEndIsSet(trans,TRUE);
  }

  return trans->end;
}

int Transcript_setCodingRegionEnd(Transcript *trans, int end) {
  trans->codingRegionEnd = end;
  Transcript_setCodingRegionEndIsSet(trans,TRUE);
  return trans->codingRegionEnd;
}

int Transcript_getCodingRegionEnd(Transcript *trans) {
  Translation *translation = Transcript_getTranslation(trans);

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
    qsort(Transcript_getExons(trans),Transcript_getExonCount(trans),
          sizeof(void *),Exon_forwardStrandCompFunc);
  } else if (strand == -1) {
    qsort(Transcript_getExons(trans),Transcript_getExonCount(trans),
          sizeof(void *),Exon_reverseStrandCompFunc);
  }
}

Exon *Transcript_getStartExon(Transcript *trans) {
  return Transcript_getExonAt(trans,0); 
}

Exon *Transcript_getEndExon(Transcript *trans) {
  return Transcript_getExonAt(trans,Transcript_getExonCount(trans)-1); 
}
