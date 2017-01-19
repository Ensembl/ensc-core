/*
 * Copyright [1999-2017] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *      http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "MapperRangeSet.h"
#include <stdio.h>
#include <stdlib.h>

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

  for (up=0; up<mrs->nRange/2; up++, down--) {
    MapperRange *tmp = mrs->ranges[down];
    mrs->ranges[down] = mrs->ranges[up];
    mrs->ranges[up]   = tmp;
  }
}

void MapperRangeSet_removeGaps(MapperRangeSet *mrs) {
  int i;
  int j;
  for (i=0; i<mrs->nRange; i++) {
    MapperRange *tmp = mrs->ranges[i];
    if (tmp->rangeType == MAPPERRANGE_GAP) {
      for (j=i+1; j<mrs->nRange; j++) {
        mrs->ranges[j-1] = mrs->ranges[j];
      }
      mrs->nRange--;
      i--;
      free(tmp);
    }
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
