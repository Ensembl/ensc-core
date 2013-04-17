#ifndef __MAPPERCOORDINATE_H__
#define __MAPPERCOORDINATE_H__

#include "EnsC.h"
#include "MapperRange.h"

typedef struct MapperCoordinateStruct MapperCoordinate;

struct MapperCoordinateStruct {
  MAPPERRANGE_DATA
  IDType id;
  signed char strand;
};

MapperCoordinate *MapperCoordinate_new(IDType id, long start, long end, int strand, CoordSystem *cs, int rank);
#endif
