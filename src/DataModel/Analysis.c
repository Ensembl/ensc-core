#include "Analysis.h"

Analysis *Analysis_new() {
  Analysis *anal;

  if ((anal = (Analysis *)calloc(1,sizeof(Analysis))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for anal\n");
    return NULL;
  }

  return anal;
}
