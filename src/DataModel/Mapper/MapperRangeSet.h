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
#define MapperRangeSet_getNumRange(mrs) (mrs)->nRange
void MapperRangeSet_removeGaps(MapperRangeSet *mrs);
void MapperRangeSet_addRange(MapperRangeSet *mrs, MapperRange *range);
void MapperRangeSet_reverse(MapperRangeSet *mrs);
void MapperRangeSet_free(MapperRangeSet *mrs);


#endif
