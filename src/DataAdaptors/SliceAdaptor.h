#ifndef __SLICEADAPTOR_H__
#define __SLICEADAPTOR_H__

#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "Slice.h"

struct SliceAdaptorStruct {
  BASEADAPTOR_DATA
};

SliceAdaptor *SliceAdaptor_new(DBAdaptor *dba);
Slice *SliceAdaptor_sliceFromRow(SliceAdaptor *sa, MYSQL_ROW row);


#endif
