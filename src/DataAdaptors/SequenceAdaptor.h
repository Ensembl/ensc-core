#ifndef __SEQUENCEADAPTOR_H__
#define __SEQUENCEADAPTOR_H__

#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "RawContig.h"
#include "Sequence.h"

struct SequenceAdaptorStruct {
  BASEADAPTOR_DATA
};

SequenceAdaptor *SequenceAdaptor_new(DBAdaptor *dba);
char *SequenceAdaptor_fetchByRawContigStartEndStrand(SequenceAdaptor *sa,
                                                     RawContig *rc,
                                                     int start,
                                                     int end,
                                                     char strand);
char *SequenceAdaptor_fetchByAssemblyLocation(SequenceAdaptor *sa,
          int chrStart, int chrEnd, int strand, long chrId, char *assemblyType);
char *SequenceAdaptor_fetchBySliceStartEndStrand(SequenceAdaptor *sa,
                                                 Slice *slice, int start, int end,
                                                 int strand);

          

#endif
