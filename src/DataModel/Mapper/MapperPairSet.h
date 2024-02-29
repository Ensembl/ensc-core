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

#ifndef __MAPPERPAIRSET_H__
#define __MAPPERPAIRSET_H__

#include "MapperPair.h"
#include "Vector.h"

typedef struct MapperPairSetStruct MapperPairSet;

struct MapperPairSetStruct {
  MapperPair **pairs;
  int nPair;
};

MapperPairSet *MapperPairSet_new();

#define MapperPairSet_getPairAt(mps,ind) (mps)->pairs[(ind)]

#define MapperPairSet_getNumPair(mps) (mps)->nPair

void MapperPairSet_sort(MapperPairSet *set, int sInd);
void MapperPairSet_free(MapperPairSet *mps, int freePairsFlag);

void MapperPairSet_addPair(MapperPairSet *mps, MapperPair *pair);

MapperPair *MapperPairSet_removePairAt(MapperPairSet *mps, int ind);

Vector *MapperPairSet_getFromIds(MapperPairSet *mps);
Vector *MapperPairSet_getToIds(MapperPairSet *mps);


#endif
