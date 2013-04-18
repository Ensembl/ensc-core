#ifndef __MAPPERPAIRSET_H__
#define __MAPPERPAIRSET_H__

#include "MapperPair.h"
#include "Vector.h"

typedef struct MapperPairSetStruct MapperPairSet;

struct MapperPairSetStruct {
  MapperPair **pairs;
  int nPair;
};

MapperPairSet *MapperPairSet_new();

#define MapperPairSet_getPairAt(mps,ind) (mps)->pairs[(ind)]

#define MapperPairSet_getNumPair(mps) (mps)->nPair

void MapperPairSet_sort(MapperPairSet *set, int sInd);
void MapperPairSet_free(MapperPairSet *mps, int freePairsFlag);

void MapperPairSet_addPair(MapperPairSet *mps, MapperPair *pair);

MapperPair *MapperPairSet_removePairAt(MapperPairSet *mps, int ind);

Vector *MapperPairSet_getFromIds(MapperPairSet *mps);
Vector *MapperPairSet_getToIds(MapperPairSet *mps);


#endif
