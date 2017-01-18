/*
 * Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
 * Copyright [2016-2017] EMBL-European Bioinformatics Institute
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

#ifndef __PREDICTIONEXONADAPTOR_H__
#define __PREDICTIONEXONADAPTOR_H__

#include "BaseFeatureAdaptor.h"
#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "PredictionExon.h"
#include "PredictionTranscript.h"

struct PredictionExonAdaptorStruct {
  BASEFEATUREADAPTOR_DATA
};

PredictionExonAdaptor *PredictionExonAdaptor_new(DBAdaptor *dba);
NameTableType *PredictionExonAdaptor_getTables();
char **PredictionExonAdaptor_getColumns();
char *PredictionExonAdaptor_finalClause();
Vector *PredictionExonAdaptor_fetchAllByPredictionTranscript(PredictionExonAdaptor *pea, PredictionTranscript *transcript);
Vector *PredictionExonAdaptor_listDbIDs(PredictionExonAdaptor *pea, int ordered);
Vector *PredictionExonAdaptor_objectsFromStatementHandle(PredictionExonAdaptor *pea, StatementHandle *sth, AssemblyMapper *assMapper, Slice *destSlice);
IDType PredictionExonAdaptor_store(PredictionExonAdaptor *pea, PredictionExon *pExon, IDType ptId, int rank);


#define PredictionExonAdaptor_genericFetch(pea, constraint, mapper, slice) \
      BaseAdaptor_genericFetch((BaseAdaptor *)(pea), (constraint), (mapper), (slice))

#define PredictionExonAdaptor_fetchAllBySlice(pea,slice)  \
   BaseFeatureAdaptor_fetchAllBySliceConstraint((BaseFeatureAdaptor *)(pea), (slice), NULL, NULL)

#define PredictionExonAdaptor_fetchAllBySliceConstraint(pea,slice,constraint,logicName)  \
   BaseFeatureAdaptor_fetchAllBySliceConstraint((BaseFeatureAdaptor *)(pea), (slice), (constraint), (logicName))


#endif
