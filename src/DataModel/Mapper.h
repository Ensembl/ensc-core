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
  CoordSystem toSystem;
  CoordSystem fromSystem;
  int isSorted;
};

#define Mapper_setTo(m, t) (m)->toSystem = (t)
#define Mapper_getTo(m) (m)->toSystem

#define Mapper_setFrom(m, f) (m)->fromSystem = (f)
#define Mapper_getFrom(m) (m)->fromSystem

#define Mapper_setPairHash(m, ind, h) (m)->hashes[(ind)] = (h)
#define Mapper_getPairHash(m, ind) (m)->hashes[(ind)]

#define Mapper_setIsSorted(m, i) (m)->isSorted = (i)
#define Mapper_getIsSorted(m) (m)->isSorted

Mapper *Mapper_new(CoordSystem from, CoordSystem to);

MapperRangeSet *Mapper_mapCoordinates(Mapper *m, long id, int start, int end, 
                                      int strand, CoordSystem type);

int Mapper_fastMap(Mapper *m, long id, int start, int end, int strand, 
                   CoordSystem type, MapperCoordinate *retRange);

void Mapper_addMapCoordinates(Mapper *m, long contigId, int contigStart, int contigEnd,
                              int contigOri, long chrId, int chrStart, int chrEnd);

MapperPairSet *Mapper_listPairs(Mapper *m, long id, int start, int end, CoordSystem type);

void Mapper_dump(Mapper *m, FILE *fp);

void Mapper_sort(Mapper *m);


#endif
