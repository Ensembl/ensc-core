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

AnalysisAdaptor *AnalysisAdaptor_new(DBAdaptor *dba);
Analysis *AnalysisAdaptor_fetchByDbID(AnalysisAdaptor *aa, long dbID);
Analysis *AnalysisAdaptor_fetchByLogicName(AnalysisAdaptor *aa, char *logicName);
Analysis *AnalysisAdaptor_analysisFromRow(AnalysisAdaptor *aa, MYSQL_ROW row);
Analysis **AnalysisAdaptor_fetchAll(AnalysisAdaptor *aa);





#endif
