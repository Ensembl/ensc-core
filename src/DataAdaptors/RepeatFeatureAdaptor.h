#ifndef __REPEATFEATUREADAPTOR_H__
#define __REPEATFEATUREADAPTOR_H__

#include "BaseFeatureAdaptor.h"
#include "AdaptorTypes.h"

struct RepeatFeatureAdaptorStruct {
  BASEFEATUREADAPTOR_DATA
};

RepeatFeatureAdaptor *RepeatFeatureAdaptor_new(DBAdaptor *dba);

#endif
