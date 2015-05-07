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

#ifndef __INDELCOORDINATE_H__
#define __INDELCOORDINATE_H__

#include "EnsC.h"
#include "MapperRange.h"
#include "MapperCoordinate.h"
#include "MapperGap.h"
#include "CoordSystem.h"

typedef struct IndelCoordinateStruct IndelCoordinate;

struct IndelCoordinateStruct {
  MAPPERRANGE_DATA
  int rank;
  IDType id;
  signed char strand;
  CoordSystem *coordSystem;
  long gapStart;
  long gapEnd;
};

IndelCoordinate *IndelCoordinate_new(MapperGap *gap, MapperCoordinate *coordinate);
#endif
