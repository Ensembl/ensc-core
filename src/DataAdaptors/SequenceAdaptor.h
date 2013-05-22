#ifndef __SEQUENCEADAPTOR_H__
#define __SEQUENCEADAPTOR_H__

#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "RawContig.h"

#include "Cache.h"
#include "IDHash.h"
#include "StringHash.h"
#include "LRUCache.h"

struct SequenceAdaptorStruct {
  BASEADAPTOR_DATA

  LRUCache *seqCache;
//  StringHash *seqCache;
  IDHash *rnaEditsCache;
};

SequenceAdaptor *SequenceAdaptor_new(DBAdaptor *dba);
StatementHandle *SequenceAdaptor_prepare(BaseAdaptor *ba, char *qStr, size_t len);

char *SequenceAdaptor_fetchByRawContigStartEndStrand(SequenceAdaptor *sa,
                                                     //RawContig *rc,
                                                     IDType rcId, 
                                                     int start,
                                                     int end,
                                                     char strand);
char *SequenceAdaptor_fetchByAssemblyLocation(SequenceAdaptor *sa,
          int chrStart, int chrEnd, int strand, char * chrName, char *assemblyType);
char *SequenceAdaptor_fetchBySliceStartEndStrand(SequenceAdaptor *sa,
                                                 Slice *slice, long start, long end,
                                                 int strand);

char *SequenceAdaptor_fetchBySliceStartEndStrandRecursive(SequenceAdaptor *sa,
                                                          Slice *slice, long start, long end,
                                                          int strand, int *recLev);

          

#endif
