#ifndef __PREDICTIONTRANSCRIPTADAPTOR_H__
#define __PREDICTIONTRANSCRIPTADAPTOR_H__

#include "BaseFeatureAdaptor.h"
#include "AdaptorTypes.h"

struct PredictionTranscriptAdaptorStruct {
  BASEFEATUREADAPTOR_DATA
};

PredictionTranscriptAdaptor *PredictionTranscriptAdaptor_new(DBAdaptor *dba);

#endif
