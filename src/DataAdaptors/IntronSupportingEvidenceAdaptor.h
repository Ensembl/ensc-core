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

#ifndef __INTRONSUPPORTINGEVIDENCEADAPTOR_H__
#define __INTRONSUPPORTINGEVIDENCEADAPTOR_H__

#include "BaseFeatureAdaptor.h"
#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "IntronSupportingEvidence.h"
#include "Intron.h"

struct IntronSupportingEvidenceAdaptorStruct {
  BASEFEATUREADAPTOR_DATA
};

IntronSupportingEvidenceAdaptor *IntronSupportingEvidenceAdaptor_new(DBAdaptor *dba);
NameTableType *IntronSupportingEvidenceAdaptor_getTables();
char **IntronSupportingEvidenceAdaptor_getColumns();
Vector *IntronSupportingEvidenceAdaptor_listLinkedTranscriptIds(IntronSupportingEvidenceAdaptor *isea, IntronSupportingEvidence *ise);
Vector *IntronSupportingEvidenceAdaptor_fetchAllByTranscript(IntronSupportingEvidenceAdaptor *isea, Transcript *transcript);
IDType *IntronSupportingEvidenceAdaptor_fetchFlankingExonIds(IntronSupportingEvidenceAdaptor *isea, IntronSupportingEvidence *ise, Transcript *transcript, IDType *flanks);
IDType IntronSupportingEvidenceAdaptor_store(IntronSupportingEvidenceAdaptor *isea, IntronSupportingEvidence *ise);
Vector *IntronSupportingEvidenceAdaptor_objectsFromStatementHandle(IntronSupportingEvidenceAdaptor *isea, StatementHandle *sth, AssemblyMapper *assMapper, Slice *destSlice);

void IntronSupportingEvidenceAdaptor_storeTranscriptLinkage(IntronSupportingEvidenceAdaptor *isea, IntronSupportingEvidence *sf, Transcript *transcript, IDType transcriptId); 

#define IntronSupportingEvidenceAdaptor_fetchAllByDbIDList(isea,id,slice)  \
   BaseAdaptor_fetchAllByDbIDList((BaseAdaptor *)(isea), (id), (slice))


#endif
