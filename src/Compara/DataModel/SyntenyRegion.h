#ifndef __SYNTENYREGION_H__
#define __SYNTENYREGION_H__

#include "ComparaDataModelTypes.h"
#include "EnsRoot.h"

struct SyntenyRegionStruct {
  ENSROOT_DATA
  Storable st;
  int start;
  int end;
  IDType clusterId;
  IDType dnaFragId;
};

#endif
