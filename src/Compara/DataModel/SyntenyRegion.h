#ifndef __SYNTENYREGION_H__
#define __SYNTENYREGION_H__

#include "ComparaDataModelTypes.h"
#include "EnsRoot.h"

#define FUNCSTRUCTTYPE NoTypeFuncs
struct SyntenyRegionStruct {
  ENSROOT_DATA
  Storable st;
  int start;
  int end;
  IDType clusterId;
  IDType dnaFragId;
};
#undef FUNCSTRUCTTYPE

#endif
