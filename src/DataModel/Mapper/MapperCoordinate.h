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

MapperCoordinate *MapperCoordinate_new(IDType id, int start, int end, int strand);
#endif
