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

#include "MapperPair.h"
#include <stdio.h>
#include <stdlib.h>


MapperPair *MapperPair_new() {
  MapperPair *mp;

  if ((mp = (MapperPair *)calloc(1,sizeof(MapperPair))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for mp\n");
    exit(1);
  }

  return mp;
}

void MapperPair_free(MapperPair *mp) {
  free(mp->units[0]);
  free(mp->units[1]);
  free(mp);
}
