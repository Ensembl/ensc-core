#ifndef __HOMOLOGY_H__
#define __HOMOLOGY_H__

#include "ComparaDataModelTypes.h"
#include "EnsRoot.h"

#define FUNCSTRUCTTYPE NoTypeFuncs
struct HomologyStruct {
  ENSROOT_DATA
  char *species;
  char *stableId;
  int chrStart;
  int chrEnd;
  char *chrName;
};
#undef FUNCSTRUCTTYPE

Homology *Homology_new();

#endif
