#ifndef __SIMPLEFEATUREADAPTOR_H__
#define __SIMPLEFEATUREADAPTOR_H__

#include "BaseFeatureAdaptor.h"
#include "AdaptorTypes.h"
#include "SimpleFeature.h"

struct SimpleFeatureAdaptorStruct {
  BASEFEATUREADAPTOR_DATA
};

SimpleFeatureAdaptor *SimpleFeatureAdaptor_new(DBAdaptor *dba);

int SimpleFeatureAdaptor_store(BaseFeatureAdaptor *baf, Vector *features);
NameTableType *SimpleFeatureAdaptor_getTables(void);
char **SimpleFeatureAdaptor_getColumns(void);
Vector *SimpleFeatureAdaptor_objectsFromStatementHandle(BaseFeatureAdaptor *baf,
                                                     StatementHandle *sth,
                                                     AssemblyMapper *mapper,
                                                     Slice *slice);



#define SimpleFeatureAdaptor_fetchByDbID(sfa, id) BaseFeatureAdaptor_fetchByDbID((BaseFeatureAdaptor *)(sfa), (id))
#define SimpleFeatureAdaptor_fetchAllBySliceAndScore(sfa, slice, score, lname) \
          BaseFeatureAdaptor_fetchAllBySliceAndScore((BaseFeatureAdaptor *)(sfa), (slice), (score), (lname))
#define SimpleFeatureAdaptor_fetchAllByRawContigAndScore(sfa, rc, score, lname) \
          BaseFeatureAdaptor_fetchAllByRawContigAndScore((BaseFeatureAdaptor *)(sfa), (rc), (score), (lname))

#endif
