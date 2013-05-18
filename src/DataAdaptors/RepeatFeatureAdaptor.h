#ifndef __REPEATFEATUREADAPTOR_H__
#define __REPEATFEATUREADAPTOR_H__

#include "BaseFeatureAdaptor.h"
#include "AdaptorTypes.h"

struct RepeatFeatureAdaptorStruct {
  BASEFEATUREADAPTOR_DATA
};

RepeatFeatureAdaptor *RepeatFeatureAdaptor_new(DBAdaptor *dba);

int RepeatFeatureAdaptor_store(BaseFeatureAdaptor *bfa, Vector *features);
char *RepeatFeatureAdaptor_defaultWhereClause();
NameTableType *RepeatFeatureAdaptor_getTables();
char **RepeatFeatureAdaptor_getColumns();
Vector *RepeatFeatureAdaptor_fetchAllBySlice(RepeatFeatureAdaptor *rfa, Slice *slice, char *logicName, Vector *repeatTypes);
Vector *RepeatFeatureAdaptor_objectsFromStatementHandle(BaseFeatureAdaptor *bfa,
                                                     StatementHandle *sth,
                                                     AssemblyMapper *mapper,
                                                     Slice *slice);


#define RepeatFeatureAdaptor_fetchByDbID(rfa, id) BaseFeatureAdaptor_fetchByDbID((BaseFeatureAdaptor *)(rfa), (id))
#define RepeatFeatureAdaptor_fetchAllBySliceConstraint(rfa, slice, constraint, lname) \
          BaseFeatureAdaptor_fetchAllBySliceConstraint((BaseFeatureAdaptor *)(rfa), (slice), (constraint), (lname))
#define RepeatFeatureAdaptor_fetchAllByRawContig(rfa, contig, lname) \
          BaseFeatureAdaptor_fetchAllByRawContig((BaseFeatureAdaptor *)(rfa), (contig), (lname))

#endif
