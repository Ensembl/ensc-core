#include "SeqFeature.h"

char *SeqFeature_setSeqName(SeqFeature *sf, char *seqName) {
  if ((sf->seqName = (char *)malloc(strlen(seqName))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for seqfeature seqName\n");
    return NULL;
  }

  strcpy(sf->seqName,seqName);

  return sf->seqName;
}
