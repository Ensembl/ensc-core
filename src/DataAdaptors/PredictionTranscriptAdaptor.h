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

#ifndef __PREDICTIONTRANSCRIPTADAPTOR_H__
#define __PREDICTIONTRANSCRIPTADAPTOR_H__

#include "BaseFeatureAdaptor.h"
#include "AdaptorTypes.h"
#include "PredictionTranscript.h"

struct PredictionTranscriptAdaptorStruct {
  BASEFEATUREADAPTOR_DATA
};

PredictionTranscriptAdaptor *PredictionTranscriptAdaptor_new(DBAdaptor *dba);
NameTableType *PredictionTranscriptAdaptor_getTables(void);
char **PredictionTranscriptAdaptor_getColumns(void);
Vector *PredictionTranscriptAdaptor_fetchAllBySlice(PredictionTranscriptAdaptor *pta, Slice *slice, char *logicName, int loadExons);
Vector *PredictionTranscriptAdaptor_objectsFromStatementHandle(PredictionTranscriptAdaptor *pta,
                                                               StatementHandle *sth,
                                                               AssemblyMapper *assMapper,
                                                               Slice *destSlice);
int PredictionTranscriptAdaptor_store(PredictionTranscriptAdaptor *pta, Vector *preTranscripts);


#define PredictionTranscriptAdaptor_genericFetch(pta, constraint, mapper, slice) \
      BaseAdaptor_genericFetch((BaseAdaptor *)(pta), (constraint), (mapper), (slice))

#define PredictionTranscriptAdaptor_fetchByDbID(pta, id) BaseFeatureAdaptor_fetchByDbID((BaseFeatureAdaptor *)(pta), (id))


#endif
