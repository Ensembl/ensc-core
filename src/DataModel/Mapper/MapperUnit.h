#ifndef __MAPPERUNIT_H__
#define __MAPPERUNIT_H__

#include "EnsC.h"

typedef struct MapperUnitStruct MapperUnit;

struct MapperUnitStruct {
  int start;
  int end;
  int ori;
  IDType id;
};

MapperUnit *MapperUnit_new();

#endif
