#ifndef __MAPPERRANGESET_H__
#define __MAPPERRANGESET_H__

#include "MapperRange.h"

typedef struct MapperRangeSetStruct MapperRangeSet;

struct MapperRangeSetStruct {
  MapperRange **ranges;
  int nRange;
};

MapperRangeSet *MapperRangeSet_new();
#define MapperRangeSet_getRangeAt(mrs,ind) (mrs)->ranges[(ind)]

#endif
