#ifndef __REPEATCONSENSUS_H__
#define __REPEATCONSENSUS_H__

#include "Storable.h"
#include "EcoString.h"
#include "EnsRoot.h"

OBJECTFUNC_TYPES(RepeatConsensus)

typedef struct RepeatConsensusFuncsStruct {
  OBJECTFUNCS_DATA(RepeatConsensus)
} RepeatConsensusFuncs;


#define FUNCSTRUCTTYPE RepeatConsensusFuncs
struct RepeatConsensusStruct {
  ENSROOT_DATA
  Storable st;
  ECOSTRING repeatClass;
  ECOSTRING consensus;
  ECOSTRING name;
  ECOSTRING repeatType;
  int length;
};
#undef FUNCSTRUCTTYPE

#define RepeatConsensus_setDbID(rcs,id) Storable_setDbID(&((rcs)->st),(id))
#define RepeatConsensus_getDbID(rcs) Storable_getDbID(&((rcs)->st))

#define RepeatConsensus_setAdaptor(rcs,ad) Storable_setAdaptor(&((rcs)->st),(ad))
#define RepeatConsensus_getAdaptor(rcs) Storable_getAdaptor(&((rcs)->st))

ECOSTRING RepeatConsensus_setRepeatClass(RepeatConsensus *rc, char *class);
#define RepeatConsensus_getRepeatClass(rcs) (rcs)->repeatClass

ECOSTRING RepeatConsensus_setRepeatType(RepeatConsensus *rc, char *type);
#define RepeatConsensus_getRepeatType(rcs) (rcs)->repeatType

ECOSTRING RepeatConsensus_setConsensus(RepeatConsensus *rc, char *cons);
#define RepeatConsensus_getConsensus(rcs) (rcs)->consensus

ECOSTRING RepeatConsensus_setName(RepeatConsensus *rc, char *name);
#define RepeatConsensus_getName(rcs) (rcs)->name

#define RepeatConsensus_setLength(rcs,len) (rcs)->length = (len)
#define RepeatConsensus_getLength(rcs) (rcs)->length

RepeatConsensus *RepeatConsensus_new();

void RepeatConsensus_free(RepeatConsensus *rc);

#ifdef __REPEATCONSENSUS_MAIN__
  RepeatConsensusFuncs
    repeatConsensusFuncs = {
                    RepeatConsensus_free,
                    NULL, // shallowCopy
                    NULL // deepCopy
                   };
#else
  extern RepeatConsensusFuncs repeatConsensusFuncs;
#endif


#endif
