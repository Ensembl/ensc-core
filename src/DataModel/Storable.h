#ifndef __STORABLE_H__
#define __STORABLE_H__

#include "BaseAdaptor.h"

typedef struct StorableStruct {
  long     dbID;
  BaseAdaptor *adaptor;
} Storable;

#define Storable_getDbID(st) (st)->dbID
#define Storable_setDbID(st,id) (st)->dbID = id
#define Storable_getAdaptor(st) (st)->adaptor
#define Storable_setAdaptor(st,ad) (st)->adaptor = ad

#endif
