#ifndef __CLONEADAPTOR_H__
#define __CLONEADAPTOR_H__

#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "Clone.h"

struct CloneAdaptorStruct {
  BASEADAPTOR_DATA
};

CloneAdaptor *CloneAdaptor_new(DBAdaptor *dba);
Clone *CloneAdaptor_cloneFromRow(CloneAdaptor *sa, ResultRow *row);


#endif
