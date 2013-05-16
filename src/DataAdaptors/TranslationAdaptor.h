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
//int TranslationAdaptor_getStableEntryInfo(TranslationAdaptor *ta, Translation *translation);
//Translation *TranslationAdaptor_fetchByDbID(TranslationAdaptor *ta, IDType dbID, Transcript *transcript);
//IDType TranslationAdaptor_store(TranslationAdaptor *ta, Translation *translation);

Vector *TranslationAdaptor_fetchAllAlternativeByTranscript(TranslationAdaptor *tlna, Transcript *transcript);
Translation *TranslationAdaptor_fetchByTranscript(TranslationAdaptor *tlna, Transcript *transcript);
Vector *TranslationAdaptor_listDbIDs(TranslationAdaptor *tlna, int ordered);
Vector *TranslationAdaptor_listStableIDs(TranslationAdaptor *tlna);
Vector *TranslationAdaptor_listDbIDsPriv(TranslationAdaptor *tlna, char *table, char *column, int ordered);
Translation *TranslationAdaptor_fetchByStableId(TranslationAdaptor *tlna, char *stableId);
Vector *TranslationAdaptor_fetchAllByTranscriptList(TranslationAdaptor *tlna, Vector *transcripts);
Translation *TranslationAdaptor_translationFromResultRow(TranslationAdaptor *tlna, ResultRow *row, Transcript *transcript);
Vector *TranslationAdaptor_fetchAll(TranslationAdaptor *tlna);

#define TranslationAdaptor_isMultiSpecies(ba) (0)

#endif
