#include "RepeatConsensus.h"
#include "EcoString.h"

RepeatConsensus *RepeatConsensus_new() {
  RepeatConsensus *rc;

  if ((rc = (RepeatConsensus *)calloc(1,sizeof(RepeatConsensus))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for rc\n");
    return NULL;
  }

  return rc;
}

ECOSTRING RepeatConsensus_setConsensus(RepeatConsensus *rc, char *cons) {
  EcoString_copyStr(ecoSTable,&(rc->consensus), cons, 0);
  return rc->consensus;
}

ECOSTRING RepeatConsensus_setRepeatClass(RepeatConsensus *rc, char *class) {
  EcoString_copyStr(ecoSTable,&(rc->repeatClass), class, 0);
  return rc->repeatClass;
}

ECOSTRING RepeatConsensus_setName(RepeatConsensus *rc, char *name) {
  EcoString_copyStr(ecoSTable,&(rc->name), name, 0);
  return rc->name;
}
