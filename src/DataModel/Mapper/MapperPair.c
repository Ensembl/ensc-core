#include "MapperPair.h"
#include <stdio.h>

MapperPair *MapperPair_new() {
  MapperPair *mp;

  if ((mp = (MapperPair *)calloc(1,sizeof(MapperPair))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for mp\n");
    exit(1);
  }

  return mp;
}
