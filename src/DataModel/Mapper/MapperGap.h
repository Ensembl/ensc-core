#ifndef __MAPPERGAP_H__
#define __MAPPERGAP_H__

#include "MapperRange.h"

typedef struct MapperGapStruct MapperGap;

struct MapperGapStruct {
  MAPPERRANGE_DATA
};

MapperGap *MapperGap_new(int start, int end);
#endif
