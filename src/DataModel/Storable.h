#ifndef __STORABLE_H__
#define __STORABLE_H__

#include "DataModelTypes.h"
#include "BaseAdaptor.h"

struct StorableStruct {
  long     dbID;
  BaseAdaptor *adaptor;
};

#define Storable_getDbID(st) (st)->dbID
#define Storable_setDbID(st,id) (st)->dbID = id
#define Storable_getAdaptor(st) (st)->adaptor
#define Storable_setAdaptor(st,ad) (st)->adaptor = ad

#endif
