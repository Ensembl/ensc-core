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

#ifndef __GENEADAPTOR_H__
#define __GENEADAPTOR_H__

#include "BaseAdaptor.h"
#include "BaseFeatureAdaptor.h"
#include "AdaptorTypes.h"
#include "Gene.h"
#include "Slice.h"
#include "Vector.h"

#include "StringHash.h"

struct GeneAdaptorStruct {
  BASEFEATUREADAPTOR_DATA
};

GeneAdaptor *GeneAdaptor_new(DBAdaptor *dba);
NameTableType *GeneAdaptor_getTables();
char **GeneAdaptor_getColumns();
IDType GeneAdaptor_store(GeneAdaptor *ga, Gene *gene, int ignoreRelease);
NameTableType *GeneAdaptor_leftJoin();
Vector *GeneAdaptor_listDbIDs(GeneAdaptor *ga, int ordered);
Vector *GeneAdaptor_listStableIDs(GeneAdaptor *ga);
Vector *GeneAdaptor_listSeqRegionIds(GeneAdaptor *ga);
Gene *GeneAdaptor_fetchByDisplayLabel(GeneAdaptor *ga, char *label);
Vector *GeneAdaptor_fetchAllByDisplayLabel(GeneAdaptor *ga, char *label);
Gene *GeneAdaptor_fetchByStableId(GeneAdaptor *ga, char *stableId);
Vector *GeneAdaptor_fetchAllByBiotype(GeneAdaptor *ga, Vector *biotypes);
void GeneAdaptor_biotypeConstraint(GeneAdaptor *ga, Vector *biotypes, char *constraint);
int GeneAdaptor_countAllByBiotype(GeneAdaptor *ga, Vector *biotypes);
Vector *GeneAdaptor_fetchAll(GeneAdaptor *ga);
Vector *GeneAdaptor_fetchAllVersionsByStableId(GeneAdaptor *ga, char *stableId);
Gene *GeneAdaptor_fetchByExonStableId(GeneAdaptor *ga, char *stableId);
int idTypeCompFunc(const void *one, const void *two);
int GeneAdaptor_countAllBySlice(GeneAdaptor *ga, Slice *slice, Vector *biotypes, char *source);
Gene *GeneAdaptor_fetchByTranscriptId(GeneAdaptor *ga, IDType transId);
Gene *GeneAdaptor_fetchByTranscriptStableId(GeneAdaptor *ga, char *transStableId);
Gene *GeneAdaptor_fetchByTranslationStableId(GeneAdaptor *ga, char *translationStableId);
Vector *GeneAdaptor_fetchAllAltAlleles(GeneAdaptor *ga, Gene *gene);
int GeneAdaptor_isRef(GeneAdaptor *ga, IDType geneId);
Vector *GeneAdaptor_objectsFromStatementHandle(GeneAdaptor *ga, StatementHandle *sth, AssemblyMapper *assMapper, Slice *destSlice);
Vector *GeneAdaptor_fetchAllBySlice(GeneAdaptor *ga, Slice *slice, char *logicName, int loadTranscripts, char *source, char *biotype);

//int GeneAdaptor_listGeneIds(GeneAdaptor *ga, IDType **geneIds);
//int GeneAdaptor_getStableEntryInfo(GeneAdaptor *ga, Gene *gene);

#define GeneAdaptor_setSpeciesId(ga,val) BaseFeatureAdaptor_setSpeciesId((ga), (val))
#define GeneAdaptor_getSpeciesId(ga) BaseFeatureAdaptor_getSpeciesId((ga))

#define GeneAdaptor_genericFetch(ga, constraint, mapper, slice) \
      BaseAdaptor_genericFetch((BaseAdaptor *)(ga), (constraint), (mapper), (slice))

#define GeneAdaptor_genericCount(ga, constraint) \
      BaseAdaptor_genericCount((BaseAdaptor *)(ga), (constraint))

#define GeneAdaptor_fetchByDbID(ga,id)  \
   BaseFeatureAdaptor_fetchByDbID((BaseFeatureAdaptor *)(ga), (id))

#define GeneAdaptor_fetchAllBySliceConstraint(ga,slice,constraint,logicName)  \
   BaseFeatureAdaptor_fetchAllBySliceConstraint((BaseFeatureAdaptor *)(ga), (slice), (constraint), (logicName))

#define GeneAdaptor_countBySliceConstraint(ga,slice,constraint,logicName)  \
   BaseFeatureAdaptor_countBySliceConstraint((BaseFeatureAdaptor *)(ga), (slice), (constraint), (logicName))

#define GeneAdaptor_fetchAllByDbIDList(ga,idList,slice)  \
   BaseAdaptor_fetchAllByDbIDList((BaseAdaptor *)(ga), (idList), (slice))




#endif
