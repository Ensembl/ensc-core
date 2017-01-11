/*
 * Copyright [1999-2017] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *      http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

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
IDType TranslationAdaptor_store(TranslationAdaptor *tlna, Translation *translation, IDType transcriptId);

#define TranslationAdaptor_isMultiSpecies(ba) (0)

#endif
