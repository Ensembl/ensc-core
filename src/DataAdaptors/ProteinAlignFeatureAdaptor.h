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

#ifndef __PROTEINALIGNFEATUREADAPTOR_H__
#define __PROTEINALIGNFEATUREADAPTOR_H__

#include "BaseFeatureAdaptor.h"
#include "AdaptorTypes.h"

struct ProteinAlignFeatureAdaptorStruct {
  BASEFEATUREADAPTOR_DATA
};

ProteinAlignFeatureAdaptor *ProteinAlignFeatureAdaptor_new(DBAdaptor *dba);
Vector *ProteinAlignFeatureAdaptor_objectsFromStatementHandle(BaseFeatureAdaptor *baf,
                                                           StatementHandle *sth,
                                                           AssemblyMapper *assMapper,
                                                           Slice *slice);
NameTableType *ProteinAlignFeatureAdaptor_leftJoin(void); 
NameTableType *ProteinAlignFeatureAdaptor_getTables(void); 
char **ProteinAlignFeatureAdaptor_getColumns(void);
int ProteinAlignFeatureAdaptor_store(BaseFeatureAdaptor *baf, Vector *features);

#define ProteinAlignFeatureAdaptor_fetchByDbID(pafa, id) BaseFeatureAdaptor_fetchByDbID((BaseFeatureAdaptor *)(pafa), (id))
#define ProteinAlignFeatureAdaptor_fetchAllBySliceAndScore(pafa, slice, score, lname) \
          BaseFeatureAdaptor_fetchAllBySliceAndScore((BaseFeatureAdaptor *)(pafa), (slice), (score), (lname))
#define ProteinAlignFeatureAdaptor_fetchAllByRawContigAndScore(pafa, contig, score, lname) \
          BaseFeatureAdaptor_fetchAllByRawContigAndScore((BaseFeatureAdaptor *)(pafa), (contig), (score), (lname))






#endif
