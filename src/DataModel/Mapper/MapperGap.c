#include "MapperGap.h"
#include <stdio.h>
#include <stdlib.h>

MapperGap *MapperGap_new(long start, long end, int rank) {
  MapperGap *mg;

  if ((mg = (MapperGap *)calloc(1,sizeof(MapperGap))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for mg\n");
    exit(1);
  }

  mg->rangeType = MAPPERRANGE_GAP;
  mg->start = start;
  mg->end = end;
  mg->rank = rank;

  return mg;
}
