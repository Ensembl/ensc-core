#include "MapperCoordinate.h"
#include <stdio.h>

MapperCoordinate *MapperCoordinate_new(long id, int start, int end, int strand) {
  MapperCoordinate *mc;

  if ((mc = (MapperCoordinate *)calloc(1,sizeof(MapperCoordinate))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for mc\n");
    exit(1);
  }

  mc->id = id;
  mc->start = start;
  mc->end = end;
  mc->strand = strand;
  return mc;
}
