#ifndef __DBENTRYADAPTOR_H__
#define __DBENTRYADAPTOR_H__

#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "DBEntry.h"

struct DBEntryAdaptorStruct {
  BASEADAPTOR_DATA
};

IDType DBEntryAdaptor_exists(DBEntryAdaptor *dbea, DBEntry *dbe);
DBEntryAdaptor *DBEntryAdaptor_new(DBAdaptor *dba);
DBEntry *DBEntryAdaptor_fetchByDbID(DBEntryAdaptor *dbea, IDType dbID);
IDType DBEntryAdaptor_store(DBEntryAdaptor *dbea, DBEntry *exObj,
                         IDType ensObject, char *ensType);

#endif
