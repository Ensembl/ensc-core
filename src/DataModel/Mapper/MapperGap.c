#include "MapperGap.h"
#include <stdio.h>

MapperGap *MapperGap_new() {
  MapperGap *mg;

  if ((mg = (MapperGap *)calloc(1,sizeof(MapperGap))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for mg\n");
    exit(1);
  }

  mg->rangeType = MAPPERRANGE_GAP;

  return mg;
}
