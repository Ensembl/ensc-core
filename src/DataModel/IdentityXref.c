#include "IdentityXref.h"
#include <stdlib.h>
#include <stdio.h>

IdentityXref *IdentityXref_new() {
  IdentityXref *idx;

  if ((idx = (IdentityXref *)calloc(1,sizeof(IdentityXref))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for idXref\n");
    return NULL;
  }

  return idx;
}

void IdentityXref_free(IdentityXref *idx) {
  free(idx);
}
