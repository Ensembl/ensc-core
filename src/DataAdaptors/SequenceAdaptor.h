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
          

#endif
