#include "MapperUnit.h"
#include <stdio.h>
#include <stdlib.h>

MapperUnit *MapperUnit_new() {
  MapperUnit *mu;

  if ((mu = (MapperUnit *)calloc(1,sizeof(MapperUnit))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for mu\n");
    exit(1);
  }

  return mu;
}
