#ifndef __REPEATCONSENSUSADAPTOR_H__
#define __REPEATCONSENSUSADAPTOR_H__

#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "RepeatConsensus.h"

struct RepeatConsensusAdaptorStruct {
  BASEADAPTOR_DATA
};

RepeatConsensusAdaptor *RepeatConsensusAdaptor_new(DBAdaptor *dba);
RepeatConsensus *RepeatConsensusAdaptor_fetchByDbID(RepeatConsensusAdaptor *ca, int64 dbID);

#endif
