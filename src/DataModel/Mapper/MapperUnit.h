#ifndef __MAPPERUNIT_H__
#define __MAPPERUNIT_H__

typedef struct MapperUnitStruct MapperUnit;

struct MapperUnitStruct {
  int start;
  int end;
  int ori;
  long id;
};

MapperUnit *MapperUnit_new();

#endif
