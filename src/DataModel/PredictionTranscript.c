#include "PredictionTranscript.h"

PredictionTranscript *PredictionTranscript_new() {
  PredictionTranscript *transcript;

  if ((transcript = (PredictionTranscript *)calloc(1,sizeof(PredictionTranscript))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for transcripe\n");
    return NULL;
  }

  return transcript;
}

char *PredictionTranscript_setType(PredictionTranscript *t, char *type) {
  if ((t->type = (char *)malloc(strlen(type)+1)) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for transcript type\n");
    return NULL;
  }

  strcpy(t->type,type);

  return t->type;
}

void PredictionTranscript_flushExons(PredictionTranscript *trans) {
// NIY  Transcript_removeAllExons(trans);
// NIY caches 
}

void PredictionTranscript_free(PredictionTranscript *trans) {
  printf("PredictionTranscript_free not implemented\n");
}