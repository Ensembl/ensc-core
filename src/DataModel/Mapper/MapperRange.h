#ifndef __MAPPERRANGE_H__
#define __MAPPERRANGE_H__

typedef enum MapperRangeTypeEnum {
  MAPPERRANGE_NONE,
  MAPPERRANGE_GAP,
  MAPPERRANGE_COORD,
  MAPPERRANGE_INDEL
} MapperRangeType;

typedef struct MapperRangeStruct MapperRange;

#define MAPPERRANGE_DATA \
   int rangeType; \
   long start; \
   long end;

struct MapperRangeStruct {
  MAPPERRANGE_DATA
};

#endif
