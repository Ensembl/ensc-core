#ifndef __PROTEINALIGNFEATUREADAPTOR_H__
#define __PROTEINALIGNFEATUREADAPTOR_H__

#include "BaseFeatureAdaptor.h"
#include "AdaptorTypes.h"

struct ProteinAlignFeatureAdaptorStruct {
  BASEFEATUREADAPTOR_DATA
};

ProteinAlignFeatureAdaptor *ProteinAlignFeatureAdaptor_new(DBAdaptor *dba);
Set *ProteinAlignFeatureAdaptor_objectsFromStatementHandle(BaseFeatureAdaptor *baf,
                                                           StatementHandle *sth,
                                                           AssemblyMapper *assMapper,
                                                           Slice *slice);
NameTableType *ProteinAlignFeatureAdaptor_getTables(void); 
char *ProteinAlignFeatureAdaptor_getColumns(void);
int ProteinAlignFeatureAdaptor_store(BaseFeatureAdaptor *baf, Set *features);

#define ProteinAlignFeatureAdaptor_fetchByDbID(pafa, id) BaseFeatureAdaptor_fetchByDbID((BaseFeatureAdaptor *)(pafa), (id))
#define ProteinAlignFeatureAdaptor_fetchAllBySliceAndScore(pafa, slice, score, lname) \
          BaseFeatureAdaptor_fetchAllBySliceAndScore((BaseFeatureAdaptor *)(pafa), (slice), (score), (lname))






#endif
