/*
 * Copyright [1999-2024] EMBL-European Bioinformatics Institute
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

#ifndef __TRANSCRIPTADAPTOR_H__
#define __TRANSCRIPTADAPTOR_H__

#include "BaseFeatureAdaptor.h"
#include "BaseFeatureAdaptor.h"
#include "AdaptorTypes.h"
#include "Transcript.h"

struct TranscriptAdaptorStruct {
  BASEFEATUREADAPTOR_DATA
};

TranscriptAdaptor *TranscriptAdaptor_new(DBAdaptor *dba);
IDType             TranscriptAdaptor_store(TranscriptAdaptor *ta, Transcript *transcript, IDType geneId, IDType analysisId);
NameTableType *    TranscriptAdaptor_getTables();
char **            TranscriptAdaptor_getColumns();
NameTableType *    TranscriptAdaptor_leftJoin();
Transcript *       TranscriptAdaptor_fetchByStableId(TranscriptAdaptor *ta, char *stableId);
Vector *           TranscriptAdaptor_fetchAll(TranscriptAdaptor *ta);
Vector *           TranscriptAdaptor_fetchAllVersionsByStableId(TranscriptAdaptor *ta, char *stableId);
Transcript *       TranscriptAdaptor_fetchByTranslationStableId(TranscriptAdaptor *ta, char *translationStableId);
Transcript *       TranscriptAdaptor_fetchByTranslationId(TranscriptAdaptor *ta, IDType translationId);
Vector *           TranscriptAdaptor_fetchAllByGene(TranscriptAdaptor *ta, Gene *gene);
Vector *           TranscriptAdaptor_fetchAllBySlice(TranscriptAdaptor *ta, Slice *slice, int loadExons, char *logicName, char *inputConstraint);
Transcript *       TranscriptAdaptor_fetchByDisplayLabel(TranscriptAdaptor *ta, char *label);
Vector *           TranscriptAdaptor_fetchByExonStableId(TranscriptAdaptor *ta, char *stableId);
Vector *           TranscriptAdaptor_fetchAllByBiotype(TranscriptAdaptor *ta, Vector *biotypes);
void               TranscriptAdaptor_biotypeConstraint(TranscriptAdaptor *ta, Vector *biotypes, char *constraint);
int                TranscriptAdaptor_isTranscriptCanonical(TranscriptAdaptor *ta, Transcript *transcript);
Vector *           TranscriptAdaptor_listDbIDs(TranscriptAdaptor *ta, int ordered);
Vector *           TranscriptAdaptor_listStableIDs(TranscriptAdaptor *ta);
Vector *           TranscriptAdaptor_objectsFromStatementHandle(TranscriptAdaptor *ta, StatementHandle *sth, AssemblyMapper *assMapper, Slice *destSlice);

#define TranscriptAdaptor_genericFetch(ta, constraint, mapper, slice) \
      BaseAdaptor_genericFetch((BaseAdaptor *)(ta), (constraint), (mapper), (slice))

#define TranscriptAdaptor_fetchByDbID(ta,id)  \
   BaseFeatureAdaptor_fetchByDbID((BaseFeatureAdaptor *)(ta), (id))

#define TranscriptAdaptor_fetchAllBySliceConstraint(ta,slice,constraint,logicName)  \
   BaseFeatureAdaptor_fetchAllBySliceConstraint((BaseFeatureAdaptor *)(ta), (slice), (constraint), (logicName))





#endif
