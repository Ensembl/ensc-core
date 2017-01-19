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

#include "IndelCoordinate.h"
#include <stdio.h>
#include <stdlib.h>

IndelCoordinate *IndelCoordinate_new(MapperGap *gap, MapperCoordinate *coordinate) {
  IndelCoordinate *ic;

  if ((ic = (IndelCoordinate *)calloc(1,sizeof(IndelCoordinate))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for ic\n");
    exit(1);
  }

  ic->id          = coordinate->id;
  ic->start       = coordinate->start;
  ic->end         = coordinate->end;
  ic->strand      = coordinate->strand;
  ic->coordSystem = coordinate->coordSystem;
  ic->gapStart    = gap->start;
  ic->gapEnd      = gap->end;
  ic->rangeType   = MAPPERRANGE_INDEL;

  return ic;
}
