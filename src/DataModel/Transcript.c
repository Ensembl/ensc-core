#include "Transcript.h"
#include "TranscriptAdaptor.h"

Transcript *Transcript_new() {
  Transcript *transcript;

  if ((transcript = (Transcript *)calloc(1,sizeof(Transcript))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for transcripe\n");
    return NULL;
  }

  return transcript;
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

  if (StableIdInfo_getVersion(&(transcript->si)) == NULL && ta) {
    TranscriptAdaptor_getStableEntryInfo(ta,transcript);
  }
  return StableIdInfo_getVersion(&(transcript->si));
}

Transcript *Transcript_transform(Transcript *trans, IDHash *exonTransforms) {
  int i;
  Set *mappedExonSet = Set_new();
  
  for (i=0;i<Transcript_getExonCount(trans);i++) {
    Exon *exon = (Exon *)Transcript_getExonAt(trans,i);

    // the old exon was successfully remapped then store the new exon
    if ( IDHash_contains(exonTransforms,(long)exon)) {
      Set_addElement(mappedExonSet,IDHash_getValue(exonTransforms, (long)exon));
    }
    // but for the case where the exon was unable to be mapped, as it
    // was outside the bounds of the slice, include the original exon.
    else {
      Set_addElement(mappedExonSet,exon);
    }
  }

#ifdef DONE
  //Flush the exons and all related internal caches
  $self->flush_Exons();

  // attach the new list of exons to the transcript
  push @{$self->{'_trans_exon_array'}},@mapped_list_of_exons;
#endif

  if ( Transcript_getTranslation(trans)) {
    Translation_transform(Transcript_getTranslation(trans), exonTransforms);
  }
}
