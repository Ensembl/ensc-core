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

#ifndef __EXONADAPTOR_H__
#define __EXONADAPTOR_H__

#include "BaseFeatureAdaptor.h"
#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "Exon.h"

struct ExonAdaptorStruct {
  BASEFEATUREADAPTOR_DATA
};

ExonAdaptor *ExonAdaptor_new(DBAdaptor *dba);
Exon *ExonAdaptor_fetchByDbID(ExonAdaptor *ea, IDType dbID);
NameTableType *ExonAdaptor_getTables();
char **ExonAdaptor_getColumns();
char *ExonAdaptor_finalClause();
Vector *ExonAdaptor_fetchAll(ExonAdaptor *ea);
Exon *ExonAdaptor_fetchByStableId(ExonAdaptor *ea, char *stableId);
Vector *ExonAdaptor_fetchAllVersionsByStableId(ExonAdaptor *ea, char *stableId);
Vector *ExonAdaptor_fetchAllByTranscript(ExonAdaptor *ea, Transcript *transcript);
IDType ExonAdaptor_store(ExonAdaptor *ea, Exon *exon);
Vector *ExonAdaptor_listDbIDs(ExonAdaptor *ea, int ordered);
Vector *ExonAdaptor_listStableIDs(ExonAdaptor *ea);
Vector *ExonAdaptor_objectsFromStatementHandle(ExonAdaptor *ea, StatementHandle *sth, AssemblyMapper *assMapper, Slice *destSlice);


#define ExonAdaptor_genericFetch(ea, constraint, mapper, slice) \
      BaseAdaptor_genericFetch((BaseAdaptor *)(ea), (constraint), (mapper), (slice))

#define ExonAdaptor_fetchByDbID(ea,id)  \
   BaseFeatureAdaptor_fetchByDbID((BaseFeatureAdaptor *)(ea), (id))

#define ExonAdaptor_fetchAllBySlice(ea,slice)  \
   BaseFeatureAdaptor_fetchAllBySliceConstraint((BaseFeatureAdaptor *)(ea), (slice), NULL, NULL)

#define ExonAdaptor_fetchAllBySliceConstraint(ea,slice,constraint,logicName)  \
   BaseFeatureAdaptor_fetchAllBySliceConstraint((BaseFeatureAdaptor *)(ea), (slice), (constraint), (logicName))


#endif
