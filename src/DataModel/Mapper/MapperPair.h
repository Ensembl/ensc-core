#ifndef __MAPPERPAIR_H__
#define __MAPPERPAIR_H__

#include "MapperDefs.h"
#include "MapperUnit.h"

typedef struct MapperPairStruct MapperPair;

struct MapperPairStruct {
  MapperUnit *units[2];
  int ori;
  int isIndel;
};

MapperPair *MapperPair_new();

#define MapperPair_setUnit(mp, ind, p) (mp)->units[(ind)] = (p)
#define MapperPair_getUnit(mp, ind) (mp)->units[(ind)]

#define MapperPair_isIndel(mp) (mp)->isIndel

void MapperPair_free(MapperPair *mp);


#endif
