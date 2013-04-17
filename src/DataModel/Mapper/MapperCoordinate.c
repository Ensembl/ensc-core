#include "MapperCoordinate.h"
#include <stdio.h>
#include <stdlib.h>

MapperCoordinate *MapperCoordinate_new(IDType id, long start, long end, int strand, CoordSystem *cs, int rank) {
  MapperCoordinate *mc;

  if ((mc = (MapperCoordinate *)calloc(1,sizeof(MapperCoordinate))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for mc\n");
    exit(1);
  }

  mc->id = id;
  mc->start = start;
  mc->end = end;
  mc->strand = strand;
  mc->coordSystem = cs;
  mc->rank = rank;
  mc->rangeType = MAPPERRANGE_COORD;
  return mc;
}
