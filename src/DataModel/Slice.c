#include "Slice.h"

Slice *Slice_new() {
  Slice *slice;

  if ((slice = (Slice *)calloc(1,sizeof(Slice))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for slice\n");
    return NULL;
  }

  return slice;
}
