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

char *Sequence_setName(Sequence *seq, char *name) {
  if ((seq->name = (char *)malloc(strlen(name)+1)) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for seq name\n");
    return NULL;
  }

  strcpy(seq->name,name);

  return seq->name;
}
