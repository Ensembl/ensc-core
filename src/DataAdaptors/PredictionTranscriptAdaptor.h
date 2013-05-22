#ifndef __PREDICTIONTRANSCRIPTADAPTOR_H__
#define __PREDICTIONTRANSCRIPTADAPTOR_H__

#include "BaseFeatureAdaptor.h"
#include "AdaptorTypes.h"
#include "PredictionTranscript.h"

struct PredictionTranscriptAdaptorStruct {
  BASEFEATUREADAPTOR_DATA
};

PredictionTranscriptAdaptor *PredictionTranscriptAdaptor_new(DBAdaptor *dba);
NameTableType *PredictionTranscriptAdaptor_getTables(void);
char **PredictionTranscriptAdaptor_getColumns(void);
Vector *PredictionTranscriptAdaptor_fetchAllBySlice(PredictionTranscriptAdaptor *pta, Slice *slice, char *logicName, int loadExons);
Vector *PredictionTranscriptAdaptor_objectsFromStatementHandle(PredictionTranscriptAdaptor *pta,
                                                               StatementHandle *sth,
                                                               AssemblyMapper *assMapper,
                                                               Slice *destSlice);
void PredictionTranscriptAdaptor_store(PredictionTranscriptAdaptor *pta, Vector *preTranscripts);


#define PredictionTranscriptAdaptor_genericFetch(pta, constraint, mapper, slice) \
      BaseAdaptor_genericFetch((BaseAdaptor *)(pta), (constraint), (mapper), (slice))

#define PredictionTranscriptAdaptor_fetchByDbID(pta, id) BaseFeatureAdaptor_fetchByDbID((BaseFeatureAdaptor *)(pta), (id))


#endif
