#ifndef __TRANSLATIONADAPTOR_H__
#define __TRANSLATIONADAPTOR_H__

#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "Translation.h"
#include "Transcript.h"

struct TranslationAdaptorStruct {
  BASEADAPTOR_DATA
};

TranslationAdaptor *TranslationAdaptor_new(DBAdaptor *dba);
int TranslationAdaptor_getStableEntryInfo(TranslationAdaptor *ta, Translation *translation);
Translation *TranslationAdaptor_fetchByDbID(TranslationAdaptor *ta, IDType dbID, Transcript *transcript);



#endif
