#ifndef __MAPPERRANGE_H__
#define __MAPPERRANGE_H__

typedef enum MapperRangeTypeEnum {
  MAPPERRANGE_GAP,
  MAPPERRANGE_COORD
} MapperRangeType;

typedef struct MapperRangeStruct MapperRange;

#define MAPPERRANGE_DATA \
   int rangeType; \
   int start; \
   int end;

struct MapperRangeStruct {
  MAPPERRANGE_DATA
};

#endif
