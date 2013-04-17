#include "CoordPair.h"
#include <stdio.h>
#include <stdlib.h>

CoordPair *CoordPair_new(long start, long end) {
  CoordPair *cp;

  if ((cp = (CoordPair *)calloc(1,sizeof(CoordPair))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for cp\n");
    exit(1);
  }

  cp->start = start;
  cp->end = end;

  return cp;
}
