#ifndef __TRANSCRIPTADAPTOR_H__
#define __TRANSCRIPTADAPTOR_H__

#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "Transcript.h"

struct TranscriptAdaptorStruct {
  BASEADAPTOR_DATA
};

TranscriptAdaptor *TranscriptAdaptor_new(DBAdaptor *dba);
int TranscriptAdaptor_getStableEntryInfo(TranscriptAdaptor *ta, Transcript *transcript);
Transcript *TranscriptAdaptor_fetchByDbID(TranscriptAdaptor *ta, int64 dbID);

#endif
