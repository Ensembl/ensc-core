#ifndef __MAPPERCOORDINATE_H__
#define __MAPPERCOORDINATE_H__

#include "EnsC.h"
#include "MapperRange.h"
#include "CoordSystem.h"

typedef struct MapperCoordinateStruct MapperCoordinate;

struct MapperCoordinateStruct {
  MAPPERRANGE_DATA
  int rank;
  IDType id;
  signed char strand;
  CoordSystem *coordSystem;
};

MapperCoordinate *MapperCoordinate_new(IDType id, long start, long end, int strand, CoordSystem *cs, int rank);
#endif
