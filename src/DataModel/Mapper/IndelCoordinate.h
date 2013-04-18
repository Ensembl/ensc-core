#ifndef __INDELCOORDINATE_H__
#define __INDELCOORDINATE_H__

#include "EnsC.h"
#include "MapperRange.h"
#include "MapperCoordinate.h"
#include "MapperGap.h"
#include "CoordSystem.h"

typedef struct IndelCoordinateStruct IndelCoordinate;

struct IndelCoordinateStruct {
  MAPPERRANGE_DATA
  int rank;
  IDType id;
  signed char strand;
  CoordSystem *coordSystem;
  long gapStart;
  long gapEnd;
};

IndelCoordinate *IndelCoordinate_new(MapperGap *gap, MapperCoordinate *coordinate);
#endif
