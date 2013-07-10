#ifndef __CACHINGSEQUENCEADAPTOR_H__
#define __CACHINGSEQUENCEADAPTOR_H__

#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "Slice.h"

#include "LRUCache.h"

struct CachingSequenceAdaptorStruct {
  BASEADAPTOR_DATA

  LRUCache *seqCache;
};

CachingSequenceAdaptor *CachingSequenceAdaptor_new(DBAdaptor *dba);
void CachingSequenceAdaptor_clearCache(CachingSequenceAdaptor *csa);

char *CachingSequenceAdaptor_fetchBySliceStartEndStrand(CachingSequenceAdaptor *csa,
                                                 Slice *slice, long start, long end,
                                                 int strand);

#endif
