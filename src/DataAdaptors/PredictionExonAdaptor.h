#ifndef __PREDICTIONEXONADAPTOR_H__
#define __PREDICTIONEXONADAPTOR_H__

#include "BaseFeatureAdaptor.h"
#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "PredictionExon.h"
#include "PredictionTranscript.h"

struct PredictionExonAdaptorStruct {
  BASEFEATUREADAPTOR_DATA
};

PredictionExonAdaptor *PredictionExonAdaptor_new(DBAdaptor *dba);
NameTableType *PredictionExonAdaptor_getTables();
char **PredictionExonAdaptor_getColumns();
char *PredictionExonAdaptor_finalClause();
Vector *PredictionExonAdaptor_fetchAllByPredictionTranscript(PredictionExonAdaptor *pea, PredictionTranscript *transcript);
Vector *PredictionExonAdaptor_listDbIDs(PredictionExonAdaptor *pea, int ordered);
Vector *PredictionExonAdaptor_objectsFromStatementHandle(PredictionExonAdaptor *pea, StatementHandle *sth, AssemblyMapper *assMapper, Slice *destSlice);
IDType PredictionExonAdaptor_store(PredictionExonAdaptor *pea, PredictionExon *pExon, IDType ptId, int rank);


#define PredictionExonAdaptor_genericFetch(pea, constraint, mapper, slice) \
      BaseAdaptor_genericFetch((BaseAdaptor *)(pea), (constraint), (mapper), (slice))

#define PredictionExonAdaptor_fetchAllBySlice(pea,slice)  \
   BaseFeatureAdaptor_fetchAllBySliceConstraint((BaseFeatureAdaptor *)(pea), (slice), NULL, NULL)

#define PredictionExonAdaptor_fetchAllBySliceConstraint(pea,slice,constraint,logicName)  \
   BaseFeatureAdaptor_fetchAllBySliceConstraint((BaseFeatureAdaptor *)(pea), (slice), (constraint), (logicName))


#endif
