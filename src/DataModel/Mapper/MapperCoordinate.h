#ifndef __MAPPERCOORDINATE_H__
#define __MAPPERCOORDINATE_H__

#include "MapperRange.h"

typedef struct MapperCoordinateStruct MapperCoordinate;

struct MapperCoordinateStruct {
  MAPPERRANGE_DATA
  long id;
  signed char strand;
};

MapperCoordinate *MapperCoordinate_new(long id, int start, int end, int strand);
#endif
