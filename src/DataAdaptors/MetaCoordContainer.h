#ifndef __METACOORDCONTAINER_H__
#define __METACOORDCONTAINER_H__

#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "IDHash.h"
#include "StringHash.h"

struct MetaCoordContainerStruct {
  BASEADAPTOR_DATA
  IDHash *maxLenCache;
  StringHash *featureCache;
};

MetaCoordContainer *MetaCoordContainer_new(DBAdaptor *dba);

#endif
