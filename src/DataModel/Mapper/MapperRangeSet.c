#include "MapperRangeSet.h"
#include <stdio.h>

MapperRangeSet *MapperRangeSet_new() {
  MapperRangeSet *mrs;

  if ((mrs = (MapperRangeSet *)calloc(1,sizeof(MapperRangeSet))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for mrs\n");
    exit(1);
  }

  return mrs;
}

void MapperRangeSet_addRange(MapperRangeSet *mrs, MapperRange *range) {
  if (mrs->nRange == 0 || !(mrs->nRange%10)) {
    if ((mrs->ranges = (MapperRange **)realloc(mrs->ranges, 
                                               (mrs->nRange+10)*sizeof(MapperRange *))) == NULL) {
      fprintf(stderr,"ERROR: Failed allocating MapperRangeSet ranges\n");
      exit(1);
    }
  }
  mrs->ranges[mrs->nRange++] = range;
}

void MapperRangeSet_reverse(MapperRangeSet *mrs) {
  int up;
  int down = mrs->nRange-1; 
  int i;

  for (up=0; up<mrs->nRange/2; up++, down--) {
    MapperRange *tmp = mrs->ranges[down];
    mrs->ranges[down] = mrs->ranges[up];
    mrs->ranges[up]   = tmp;
  }
}

void MapperRangeSet_free(MapperRangeSet *mrs) {
  int i;

  for (i=0;i<mrs->nRange;i++) {
    free(mrs->ranges[i]);
  }
  free(mrs->ranges);
  free(mrs);
}
