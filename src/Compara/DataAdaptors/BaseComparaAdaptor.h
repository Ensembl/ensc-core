#ifndef __BASECOMPARAADAPTOR_H__
#define __BASECOMPARAADAPTOR_H__

#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "ComparaDBAdaptor.h"

#define BASECOMPARAADAPTOR_DATA BASEADAPTOR_DEF(ComparaDBAdaptor)
struct BaseComparaAdaptorStruct {
  BASECOMPARAADAPTOR_DATA
};

void BaseComparaAdaptor_init(BaseComparaAdaptor *ba, ComparaDBAdaptor *dba, int adaptorType);

#endif
