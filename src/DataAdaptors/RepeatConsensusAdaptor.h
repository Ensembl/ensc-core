#ifndef __REPEATCONSENSUSADAPTOR_H__
#define __REPEATCONSENSUSADAPTOR_H__

#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "RepeatConsensus.h"
#include "Vector.h"

struct RepeatConsensusAdaptorStruct {
  BASEADAPTOR_DATA
};

RepeatConsensusAdaptor *RepeatConsensusAdaptor_new(DBAdaptor *dba);
RepeatConsensus *RepeatConsensusAdaptor_fetchByDbID(RepeatConsensusAdaptor *ca, IDType dbID);
Vector *RepeatConsensusAdaptor_genericFetch(RepeatConsensusAdaptor *rca, char *whereClause);
RepeatConsensus *RepeatConsensusAdaptor_fetchByName(RepeatConsensusAdaptor *rca, char *name);
RepeatConsensus *RepeatConsensusAdaptor_fetchByNameAndClass(RepeatConsensusAdaptor *rca, char *name, char *class);
Vector *RepeatConsensusAdaptor_fetchByClassAndSeq(RepeatConsensusAdaptor *rca, char *class, char *seq);
int RepeatConsensusAdaptor_store(RepeatConsensusAdaptor *rca, Vector *consensi);
int RepeatConsensusAdaptor_free(RepeatConsensusAdaptor *rc);


#endif
