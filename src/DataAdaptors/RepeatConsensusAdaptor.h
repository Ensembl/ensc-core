#ifndef __REPEATCONSENSUSADAPTOR_H__
#define __REPEATCONSENSUSADAPTOR_H__

#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "RepeatConsensus.h"
#include "Set.h"

struct RepeatConsensusAdaptorStruct {
  BASEADAPTOR_DATA
};

RepeatConsensusAdaptor *RepeatConsensusAdaptor_new(DBAdaptor *dba);
RepeatConsensus *RepeatConsensusAdaptor_fetchByDbID(RepeatConsensusAdaptor *ca, int64 dbID);
Set *RepeatConsensusAdaptor_genericFetch(RepeatConsensusAdaptor *rca, char *whereClause);


#endif
