#ifndef __REPEATCONSENSUS_H__
#define __REPEATCONSENSUS_H__

#include "Storable.h"

struct RepeatConsensusStruct {
  Storable st;
};

#define RepeatConsensus_getDbID(rcs) Storable_getDbID(&((rcs)->st))
#define RepeatConsensus_setDbID(rcs,id) Storable_setDbID(&((rcs)->st),dbID)
#define RepeatConsensus_getAdaptor(rcs) Storable_getAdaptor(&((rcs)->st))
#define RepeatConsensus_setAdaptor(rcs,ad) Storable_setAdaptor(&((rcs->st),ad)

#endif
