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

#ifndef __MAPPER_H__
#define __MAPPER_H__

#include <stdio.h>

#include "IDHash.h"

#include "DataModelTypes.h"

#include "MapperDefs.h"
#include "MapperPairSet.h"
#include "MapperRangeSet.h"
#include "MapperCoordinate.h"
#include "MapperRange.h"
#include "MapperGap.h"

struct MapperStruct {
  IDHash *hashes[2];
  char *from;
  char *to;
  CoordSystem *toSystem;
  CoordSystem *fromSystem;
  int isSorted;
  int pairCount;
};

#define Mapper_getTo(m) (m)->to
#define Mapper_getFrom(m) (m)->from

#define Mapper_setToCoordSystem(m, t) (m)->toSystem = (t)
#define Mapper_getToCoordSystem(m) (m)->toSystem

#define Mapper_setFromCoordSystem(m, f) (m)->fromSystem = (f)
#define Mapper_getFromCoordSystem(m) (m)->fromSystem

#define Mapper_setPairHash(m, ind, h) (m)->hashes[(ind)] = (h)
#define Mapper_getPairHash(m, ind) (m)->hashes[(ind)]

#define Mapper_setIsSorted(m, i) (m)->isSorted = (i)
#define Mapper_getIsSorted(m) (m)->isSorted
#define Mapper_isSorted(m) (m)->isSorted

#define Mapper_setPairCount(m, pc) (m)->pairCount = (pc)
#define Mapper_getPairCount(m) (m)->pairCount
#define Mapper_incPairCount(m) ((m)->pairCount++)
#define Mapper_addToPairCount(m, n) ((m)->pairCount+=n)

#define Mapper_compareType(a, b) (strcmp((a),(b)))

Mapper *Mapper_new(char *from, char *to, CoordSystem *fromCs, CoordSystem *toCs);
char *Mapper_setFrom(Mapper *m, char *from);
char *Mapper_setTo(Mapper *m, char *to);

void Mapper_flush(Mapper *m);

MapperRangeSet *Mapper_mapCoordinates(Mapper *m, IDType id, long start, long end, 
                                      int strand, char *type);

MapperRangeSet *Mapper_fastMap(Mapper *m, IDType id, long start, long end, int strand, 
                               char *type);

void Mapper_addMapCoordinates(Mapper *m, IDType contigId, int contigStart, int contigEnd,
                              int contigOri, IDType chrId, int chrStart, int chrEnd);

MapperRangeSet *Mapper_mapInsert(Mapper *m, IDType id, long start, long end, int strand, char *type, int fastmap);

MapperPairSet *Mapper_listPairs(Mapper *m, IDType id, long start, long end, char *type);

void Mapper_dump(Mapper *m, FILE *fp);

void Mapper_sort(Mapper *m);

void Mapper_mergePairs(Mapper *m);

void Mapper_free(Mapper *mapper);


#endif
