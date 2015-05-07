/*
 * Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

#ifndef __MAPPERRANGESET_H__
#define __MAPPERRANGESET_H__

#include "MapperRange.h"

typedef struct MapperRangeSetStruct MapperRangeSet;

struct MapperRangeSetStruct {
  MapperRange **ranges;
  int nRange;
};

MapperRangeSet *MapperRangeSet_new();
#define MapperRangeSet_getRangeAt(mrs,ind) (mrs)->ranges[(ind)]
#define MapperRangeSet_getNumRange(mrs) (mrs)->nRange
void MapperRangeSet_removeGaps(MapperRangeSet *mrs);
void MapperRangeSet_addRange(MapperRangeSet *mrs, MapperRange *range);
void MapperRangeSet_reverse(MapperRangeSet *mrs);
void MapperRangeSet_free(MapperRangeSet *mrs);


#endif
