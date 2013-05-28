#ifndef __TRANSCRIPTSUPPORTINGFEATUREADAPTOR_H__
#define __TRANSCRIPTSUPPORTINGFEATUREADAPTOR_H__

#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "Transcript.h"

struct TranscriptSupportingFeatureAdaptorStruct {
  BASEADAPTOR_DATA
};

TranscriptSupportingFeatureAdaptor *TranscriptSupportingFeatureAdaptor_new(DBAdaptor *dba);
Vector *TranscriptSupportingFeatureAdaptor_fetchAllByTranscript(TranscriptSupportingFeatureAdaptor *sfa, Transcript *transcript);
void TranscriptSupportingFeatureAdaptor_store(TranscriptSupportingFeatureAdaptor *tsfa, IDType transcriptDbID, Vector *alnObjs);


#endif
