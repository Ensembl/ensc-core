#include "RawContig.h"
#include "RawContigAdaptor.h"

RawContig *RawContig_new() {
  RawContig *rc;

  if ((rc = (RawContig *)calloc(1,sizeof(RawContig))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for rc\n");
    return NULL;
  }

  rc->length = rc->emblOffset = rc->cloneId = -1;

  return rc;
}

char *RawContig_setName(RawContig *rc, char *name) {
  if ((rc->name = (char *)malloc(strlen(name))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for seqfeature name\n");
    return NULL;
  }

  strcpy(rc->name,name);

  return rc->name;
}

int RawContig_getLength(RawContig *rc) {
  RawContigAdaptor *rca = (RawContigAdaptor *)RawContig_getAdaptor(rc);

  if (rc->length == -1 && rca) {
    RawContigAdaptor_fetchAttributes(rca,rc);
  }
  return rc->length;
}

int RawContig_getEMBLOffset(RawContig *rc) {
  RawContigAdaptor *rca = (RawContigAdaptor *)RawContig_getAdaptor(rc);

  if (rc->emblOffset == -1 && rca) {
    RawContigAdaptor_fetchAttributes(rca,rc);
  }
  return rc->emblOffset;
}

long RawContig_getCloneID(RawContig *rc) {
  RawContigAdaptor *rca = (RawContigAdaptor *)RawContig_getAdaptor(rc);

  if (rc->cloneId == -1 && rca) {
    RawContigAdaptor_fetchAttributes(rca,rc);
  }
  return rc->cloneId;
}

char *RawContig_getName(RawContig *rc) {
  RawContigAdaptor *rca = (RawContigAdaptor *)RawContig_getAdaptor(rc);

  if (rc->name == NULL && rca) {
    RawContigAdaptor_fetchAttributes(rca,rc);
  }
  return rc->name;
}
