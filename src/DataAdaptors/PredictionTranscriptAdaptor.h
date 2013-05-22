#ifndef __PREDICTIONTRANSCRIPTADAPTOR_H__
#define __PREDICTIONTRANSCRIPTADAPTOR_H__

#include "BaseFeatureAdaptor.h"
#include "AdaptorTypes.h"
#include "PredictionTranscript.h"

struct PredictionTranscriptAdaptorStruct {
  BASEFEATUREADAPTOR_DATA
};

PredictionTranscriptAdaptor *PredictionTranscriptAdaptor_new(DBAdaptor *dba);
int PredictionTranscriptAdaptor_store(BaseFeatureAdaptor *bfa, Vector *features);
NameTableType *PredictionTranscriptAdaptor_getTables(void);
char *PredictionTranscriptAdaptor_getColumns(void);
Vector *PredictionTranscriptAdaptor_objectsFromStatementHandle(BaseFeatureAdaptor *bfa,
                                                       StatementHandle *sth,
                                                       AssemblyMapper *mapper,
                                                       Slice *slice);

#define PredictionTranscriptAdaptor_fetchByDbID(pta, id) BaseFeatureAdaptor_fetchByDbID((BaseFeatureAdaptor *)(pta), (id))
#define PredictionTranscriptAdaptor_fetchAllBySlice(pta, slice) \
          BaseFeatureAdaptor_fetchAllBySlice((BaseFeatureAdaptor *)(pta), (slice), NULL, NULL)


#endif
