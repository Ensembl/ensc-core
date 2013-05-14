#ifndef __EXONADAPTOR_H__
#define __EXONADAPTOR_H__

#include "BaseFeatureAdaptor.h"
#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "Exon.h"

struct ExonAdaptorStruct {
  BASEFEATUREADAPTOR_DATA
};

ExonAdaptor *ExonAdaptor_new(DBAdaptor *dba);
Exon *ExonAdaptor_fetchByDbID(ExonAdaptor *ea, IDType dbID);
NameTableType *ExonAdaptor_getTables();
char **ExonAdaptor_getColumns();
char *ExonAdaptor_finalClause();
Vector *ExonAdaptor_fetchAll(ExonAdaptor *ea);
Exon *ExonAdaptor_fetchByStableId(ExonAdaptor *ea, char *stableId);
Vector *ExonAdaptor_fetchAllVersionsByStableId(ExonAdaptor *ea, char *stableId);
Vector *ExonAdaptor_fetchAllByTranscript(ExonAdaptor *ea, Transcript *transcript);
IDType ExonAdaptor_store(ExonAdaptor *ea, Exon *exon);
Vector *ExonAdaptor_listDbIDs(ExonAdaptor *ea, int ordered);
Vector *ExonAdaptor_listStableIDs(ExonAdaptor *ea);
Vector *ExonAdaptor_objectsFromStatementHandle(ExonAdaptor *ea, StatementHandle *sth, AssemblyMapper *assMapper, Slice *destSlice);


#define ExonAdaptor_genericFetch(ea, constraint, mapper, slice) \
      BaseAdaptor_genericFetch((BaseAdaptor *)(ea), (constraint), (mapper), (slice))

#define ExonAdaptor_fetchByDbID(ea,id)  \
   BaseFeatureAdaptor_fetchByDbID((BaseFeatureAdaptor *)(ea), (id))

#define ExonAdaptor_fetchAllBySliceConstraint(ea,slice,constraint,logicName)  \
   BaseFeatureAdaptor_fetchAllBySliceConstraint((BaseFeatureAdaptor *)(ea), (slice), (constraint), (logicName))


#endif
