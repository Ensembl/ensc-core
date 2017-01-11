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

#ifndef __MAPPERPAIR_H__
#define __MAPPERPAIR_H__

#include "MapperDefs.h"
#include "MapperUnit.h"

typedef struct MapperPairStruct MapperPair;

struct MapperPairStruct {
  MapperUnit *units[2];
  int ori;
  int isIndel;
};

MapperPair *MapperPair_new();

#define MapperPair_setUnit(mp, ind, p) (mp)->units[(ind)] = (p)
#define MapperPair_getUnit(mp, ind) (mp)->units[(ind)]

#define MapperPair_isIndel(mp) (mp)->isIndel

void MapperPair_free(MapperPair *mp);


#endif
