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
RepeatConsensus *RepeatConsensusAdaptor_fetchByDbID(RepeatConsensusAdaptor *ca, IDType dbID);
Set *RepeatConsensusAdaptor_genericFetch(RepeatConsensusAdaptor *rca, char *whereClause);
RepeatConsensus *RepeatConsensusAdaptor_fetchByName(RepeatConsensusAdaptor *rca, char *name);
RepeatConsensus *RepeatConsensusAdaptor_fetchByNameAndClass(RepeatConsensusAdaptor *rca, char *name, char *class);
Set *RepeatConsensusAdaptor_fetchByClassAndSeq(RepeatConsensusAdaptor *rca, char *class, char *seq);
int RepeatConsensus_free(RepeatConsensus *rc);
int RepeatConsensusAdaptor_store(RepeatConsensusAdaptor *rca, Set *consensi);

#endif
