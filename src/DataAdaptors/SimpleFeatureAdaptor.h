#ifndef __SIMPLEFEATUREADAPTOR_H__
#define __SIMPLEFEATUREADAPTOR_H__

#include "BaseFeatureAdaptor.h"
#include "AdaptorTypes.h"
#include "SimpleFeature.h"

struct SimpleFeatureAdaptorStruct {
  BASEFEATUREADAPTOR_DATA
};

SimpleFeatureAdaptor *SimpleFeatureAdaptor_new(DBAdaptor *dba);

int SimpleFeatureAdaptor_store(BaseFeatureAdaptor *baf, Set *features);
char ***SimpleFeatureAdaptor_getTables(void);
char *SimpleFeatureAdaptor_getColumns(void);
Set *SimpleFeatureAdaptor_objectsFromStatementHandle(BaseFeatureAdaptor *baf,
                                                     StatementHandle *sth,
                                                     AssemblyMapper *mapper,
                                                     Slice *slice);



#define SimpleFeatureAdaptor_fetchByDbID(sfa, id) BaseFeatureAdaptor_fetchByDbID((BaseFeatureAdaptor *)(sfa), (id))

#endif
