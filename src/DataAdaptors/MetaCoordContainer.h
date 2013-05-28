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
Vector *MetaCoordContainer_fetchAllCoordSystemsByFeatureType(MetaCoordContainer *mcc, char *origTable);
long MetaCoordContainer_fetchMaxLengthByCoordSystemFeatureType(MetaCoordContainer *mcc, CoordSystem *cs, char *table);
void MetaCoordContainer_addFeatureType(MetaCoordContainer *mcc, CoordSystem *cs, char *table, long length);


#endif
