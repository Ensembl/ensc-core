#include "IndelCoordinate.h"
#include <stdio.h>
#include <stdlib.h>

IndelCoordinate *IndelCoordinate_new(MapperGap *gap, MapperCoordinate *coordinate) {
  IndelCoordinate *ic;

  if ((ic = (IndelCoordinate *)calloc(1,sizeof(IndelCoordinate))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for ic\n");
    exit(1);
  }

  ic->id          = coordinate->id;
  ic->start       = coordinate->start;
  ic->end         = coordinate->end;
  ic->strand      = coordinate->strand;
  ic->coordSystem = coordinate->coordSystem;
  ic->gapStart    = gap->start;
  ic->gapEnd      = gap->end;
  ic->rangeType   = MAPPERRANGE_INDEL;

  return ic;
}
