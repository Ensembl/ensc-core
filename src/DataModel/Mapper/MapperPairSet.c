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

void MapperPairSet_addPair(MapperPairSet *mps, MapperPair *pair) {
  if (mps->nPair == 0 || !(mps->nPair%10)) {
    if ((mps->pairs = (MapperPair **)realloc(mps->pairs, 
                                               (mps->nPair+10)*sizeof(MapperPair *))) == NULL) {
      fprintf(stderr,"ERROR: Failed allocating MapperPairSet pairs\n");
      exit(1);
    }
  } 
  mps->pairs[mps->nPair++] = pair;
}

static int sortInd;

int MapperPairCompFunc(const void *a, const void *b) {
  MapperPair **mp1 = (MapperPair **)a;
  MapperPair **mp2 = (MapperPair **)b;

  printf("Comparing %d %d\n",MapperPair_getUnit(*mp1,sortInd)->start,
         MapperPair_getUnit(*mp2,sortInd)->start);
  if (MapperPair_getUnit(*mp1,sortInd)->start >
      MapperPair_getUnit(*mp2,sortInd)->start) {
    return 1;
  } else if (MapperPair_getUnit(*mp1,sortInd)->start <
             MapperPair_getUnit(*mp2,sortInd)->start) {
    return -1;
  } else {
    return 0;
  }
}

void MapperPairSet_sort(MapperPairSet *mps, int sInd) {
  sortInd = sInd;
  qsort(mps->pairs,(size_t)mps->nPair,sizeof(MapperPair *),MapperPairCompFunc);
}

void MapperPairSet_free(MapperPairSet *mps, int freePairsFlag) {
  int i;
 
  if (freePairsFlag) {
    for (i=0;i<mps->nPair;i++) {
      MapperPair_free(mps->pairs[i]);
    }
  }
  free(mps->pairs);
  free(mps);
}
