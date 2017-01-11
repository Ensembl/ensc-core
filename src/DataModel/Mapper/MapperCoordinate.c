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

#include "MapperCoordinate.h"
#include <stdio.h>
#include <stdlib.h>

MapperCoordinate *MapperCoordinate_new(IDType id, long start, long end, int strand, CoordSystem *cs, int rank) {
  MapperCoordinate *mc;

  if ((mc = (MapperCoordinate *)calloc(1,sizeof(MapperCoordinate))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for mc\n");
    exit(1);
  }

  mc->id = id;
  mc->start = start;
  mc->end = end;
  mc->strand = strand;
  mc->coordSystem = cs;
  mc->rank = rank;
  mc->rangeType = MAPPERRANGE_COORD;
  return mc;
}
