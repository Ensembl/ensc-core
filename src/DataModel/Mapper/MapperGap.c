#include "MapperGap.h"
#include <stdio.h>
#include <stdlib.h>

MapperGap *MapperGap_new(int start, int end) {
  MapperGap *mg;

  if ((mg = (MapperGap *)calloc(1,sizeof(MapperGap))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for mg\n");
    exit(1);
  }

  mg->rangeType = MAPPERRANGE_GAP;
  mg->start = start;
  mg->end = end;

  return mg;
}
