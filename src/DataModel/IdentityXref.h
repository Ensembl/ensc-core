#ifndef __IDENTITYXREF_H__
#define __IDENTITYXREF_H__

#include "DataModelTypes.h"

struct IdentityXrefStruct {
  double queryIdentity;
  double targetIdentity;
};

IdentityXref *IdentityXref_new();

#define IdentityXref_setQueryIdentity(idx, qi) (idx)->queryIdentity = (qi)
#define IdentityXref_getQueryIdentity(idx) (idx)->queryIdentity

#define IdentityXref_setTargetIdentity(idx, ti) (idx)->targetIdentity = (ti)
#define IdentityXref_getTargetIdentity(idx) (idx)->targetIdentity

#endif
