#ifndef __BASEDBADAPTOR_H__
#define __BASEDBADAPTOR_H__

#include "DBConnection.h"
#include "AdaptorTypes.h"
#include "EnsC.h"

#define BASEDBADAPTOR_DATA \
  DBConnection  *dbc; \
  MetaContainer *metaContainer; \
  MetaCoordContainer *metaCoordContainer;

struct BaseDBAdaptorStruct {
  BASEDBADAPTOR_DATA
};


#define BaseDBAdaptor_prepare(dba,qStr,qLen) (dba)->dbc->prepare((dba)->dbc,(qStr),(qLen))


#endif
