#ifndef __ANALYSISADAPTOR_H__
#define __ANALYSISADAPTOR_H__

#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "Analysis.h"
#include "IDHash.h"
#include "StringHash.h"

struct AnalysisAdaptorStruct {
  BASEADAPTOR_DATA
  IDHash *analCache;
  StringHash *logicNameCache;
};

int64 AnalysisAdaptor_analysisExists(AnalysisAdaptor *aa, Analysis *anal);
AnalysisAdaptor *AnalysisAdaptor_new(DBAdaptor *dba);
Analysis *AnalysisAdaptor_fetchByDbID(AnalysisAdaptor *aa, int64 dbID);
Analysis *AnalysisAdaptor_fetchByLogicName(AnalysisAdaptor *aa, char *logicName);
Analysis *AnalysisAdaptor_analysisFromRow(AnalysisAdaptor *aa, ResultRow *row);
Analysis **AnalysisAdaptor_fetchAll(AnalysisAdaptor *aa);
int64 AnalysisAdaptor_store(AnalysisAdaptor *aa, Analysis *analysis);






#endif
