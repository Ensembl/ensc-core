#ifndef __MAPPERPAIRSET_H__
#define __MAPPERPAIRSET_H__

#include "MapperPair.h"

typedef struct MapperPairSetStruct MapperPairSet;

struct MapperPairSetStruct {
  MapperPair **pairs;
  int nPair;
};

MapperPairSet *MapperPairSet_new();

#define MapperPairSet_getPairAt(mps,ind) (mps)->pairs[(ind)]

#define MapperPairSet_getNumPair(mps) (mps)->nPair

#endif
