#define __PREDICTIONTRANSCRIPT_MAIN__
#include "PredictionTranscript.h"
#undef  __PREDICTIONTRANSCRIPT_MAIN__

#include "Exon.h"

PredictionTranscript *PredictionTranscript_new() {
  PredictionTranscript *transcript;

  if ((transcript = (PredictionTranscript *)calloc(1,sizeof(PredictionTranscript))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for transcripe\n");
    return NULL;
  }

  transcript->objectType = CLASS_PREDICTIONTRANSCRIPT;
  transcript->funcs = &predictionTranscriptFuncs;
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

int PredictionTranscript_getLength(PredictionTranscript *trans) {
    int length = 0;
    int i;

    for (i=0;i<PredictionTranscript_getExonCount(trans); i++) {
      Exon *ex = PredictionTranscript_getExonAt(trans,i);
 // Check stickies
      if (ex) length += Exon_getLength(ex);
    }
    return length;
}


