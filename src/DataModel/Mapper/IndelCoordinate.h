#ifndef __MAPPERCOORDINATE_H__
#define __MAPPERCOORDINATE_H__

#include "EnsC.h"
#include "MapperRange.h"

typedef struct IndelCoordinateStruct IndelCoordinate;

struct IndelCoordinateStruct {
  MAPPERRANGE_DATA
  IDType id;
  signed char strand;
  long gapStart;
  long gapEnd;
};

IndelCoordinate *IndelCoordinate_new(MapperGap *gap, MapperCoordinate *coordinate) {
#endif
