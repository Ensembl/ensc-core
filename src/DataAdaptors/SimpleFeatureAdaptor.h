#ifndef __SIMPLEFEATUREADAPTOR_H__
#define __SIMPLEFEATUREADAPTOR_H__

#include "BaseFeatureAdaptor.h"
#include "AdaptorTypes.h"

struct SimpleFeatureAdaptorStruct {
  BASEFEATUREADAPTOR_DATA
};

SimpleFeatureAdaptor *SimpleFeatureAdaptor_new(DBAdaptor *dba);

#endif
