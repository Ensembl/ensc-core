#include "MapperPairSet.h"
#include <stdio.h>
#include <stdlib.h>

MapperPairSet *MapperPairSet_new() {
  MapperPairSet *mps;

  if ((mps = (MapperPairSet *)calloc(1,sizeof(MapperPairSet))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for mps\n");
    exit(1);
  }

  return mps;
}

void MapperPairSet_addPair(MapperPairSet *mps, MapperPair *range) {
  if (mps->nPair == 0 || !(mps->nPair%10)) {
    if ((mps->pairs = (MapperPair **)realloc(mps->pairs, 
                                               (mps->nPair+10)*sizeof(MapperPair *))) == NULL) {
      fprintf(stderr,"ERROR: Failed allocating MapperPairSet pairs\n");
      exit(1);
    }
  } 
  mps->pairs[mps->nPair++] = range;
}


int MapperPairCompFunc(const void *a, const void *b) {
  MapperPair *mp1 = (MapperPair *)a;
  MapperPair *mp2 = (MapperPair *)b;

  return (MapperPair_getUnit(mp2,MAPPER_FROM_IND)->start - 
          MapperPair_getUnit(mp1,MAPPER_FROM_IND)->start);
}

void MapperPairSet_sort(MapperPairSet *mps) {
  qsort(mps->pairs,(size_t)mps->nPair,sizeof(MapperPair *),MapperPairCompFunc);
}
