#define __SEQUENCE_MAIN__
#include "Sequence.h"
#undef __SEQUENCE_MAIN__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

Sequence *Sequence_new() {
  Sequence *seq;

  if ((seq = (Sequence *)calloc(1,sizeof(Sequence))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for seq\n");
    return NULL;
  }

  seq->length = -1;

  seq->objectType = CLASS_SEQUENCE;

  seq->funcs = &sequenceFuncs;

  return seq;
}

ECOSTRING Sequence_setName(Sequence *seq, char *name) {
  if (EcoString_copyStr(ecoSTable,&(seq->name),name,0)) {
    fprintf(stderr,"ERROR: Failed allocating space for seq name\n");
    return NULL;
  }

  return seq->name;
}

char *Sequence_setSeq(Sequence *seq, char *seqstr) {
  if (StrUtil_copyString(&(seq->seq),seqstr,0)) {
    fprintf(stderr,"ERROR: Failed allocating space for seq\n");
    return NULL;
  }

  seq->length = strlen(seq->seq);

  return seq->seq;
}

void Sequence_freePtrs(Sequence *seq) {
  if (seq->name) free(seq->name);
// Is this OK!??
  if (seq->seq) free(seq->seq);

}
