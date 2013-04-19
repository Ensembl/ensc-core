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
void *RangeRegistry_flush(RangeRegistry *registry);
Vector *RangeRegistry_checkAndRegister(RangeRegistry *reg, IDType id, long start, long end, long rStart, long rEnd, int wantGaps);
long RangeRegistry_overlapSize(RangeRegistry *registry, IDType id, long start, long end);
Vector *RangeRegistry_getRanges(RangeRegistry *registry, IDType id);

#define RangeRegistry_getRegistry(r) (r)->registry

#endif
