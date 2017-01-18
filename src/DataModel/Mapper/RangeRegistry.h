/*
 * Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
 * Copyright [2016-2017] EMBL-European Bioinformatics Institute
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

#ifndef __RANGEREGISTRY_H__
#define __RANGEREGISTRY_H__

#include "EnsC.h"
#include "IDHash.h"
#include "Vector.h"

typedef struct RangeRegistryStruct RangeRegistry;

struct RangeRegistryStruct {
  IDHash *registry;
};

RangeRegistry *RangeRegistry_new();
void RangeRegistry_flush(RangeRegistry *registry);
Vector *RangeRegistry_checkAndRegister(RangeRegistry *reg, IDType id, long start, long end, long rStart, long rEnd, int wantGaps);
long RangeRegistry_overlapSize(RangeRegistry *registry, IDType id, long start, long end);
Vector *RangeRegistry_getRanges(RangeRegistry *registry, IDType id);

#define RangeRegistry_getRegistry(r) (r)->registry

#endif
