#define __SEQFEATURE_MAIN__
#include "SeqFeature.h"
#undef __SEQFEATURE_MAIN__

char *SeqFeature_setSeqName(SeqFeature *sf, char *seqName) {
  if ((sf->seqName = (char *)malloc(strlen(seqName)+1)) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for seqfeature seqName\n");
    return NULL;
  }

  strcpy(sf->seqName,seqName);

  return sf->seqName;
}

int SeqFeature_startCompFunc(const void *a, const void *b) {
  SeqFeature **e1 = (SeqFeature **)a;
  SeqFeature **e2 = (SeqFeature **)b;

  if (SeqFeature_getStart(*e1) > SeqFeature_getStart(*e2)) {
    return 1;
  } else if (SeqFeature_getStart(*e1) < SeqFeature_getStart(*e2)) {
    return -1;
  } else {
    return 0;
  }
}
