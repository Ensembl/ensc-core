#ifndef __HOMOLOGY_H__
#define __HOMOLOGY_H__

#include "ComparaDataModelTypes.h"
#include "EnsRoot.h"

struct HomologyStruct {
  ENSROOT_DATA
  char *species;
  char *stableId;
  int chrStart;
  int chrEnd;
  char *chrName;
};

#endif
