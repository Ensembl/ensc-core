/*
 * Copyright [1999-2024] EMBL-European Bioinformatics Institute
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

#ifndef __MAPPERCOORDINATE_H__
#define __MAPPERCOORDINATE_H__

#include "EnsC.h"
#include "MapperRange.h"
#include "CoordSystem.h"

typedef struct MapperCoordinateStruct MapperCoordinate;

struct MapperCoordinateStruct {
  MAPPERRANGE_DATA
  int rank;
  IDType id;
  signed char strand;
  CoordSystem *coordSystem;
};

MapperCoordinate *MapperCoordinate_new(IDType id, long start, long end, int strand, CoordSystem *cs, int rank);
#endif
