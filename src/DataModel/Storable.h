#ifndef __STORABLE_H__
#define __STORABLE_H__

#include "EnsC.h"
#include "DataModelTypes.h"
#include "AdaptorTypes.h"

struct StorableStruct {
  IDType     dbID;
  BaseAdaptor *adaptor;
};


#define Storable_getDbID(st) (st)->dbID
#define Storable_setDbID(st,id) (st)->dbID = id
#define Storable_getAdaptor(st) (st)->adaptor
#define Storable_setAdaptor(st,ad) (st)->adaptor = ad

int Storable_isStored(Storable *storable, DBAdaptor *db);

#endif
