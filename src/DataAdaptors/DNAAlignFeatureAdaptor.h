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

#ifndef __DNAALIGNFEATUREADAPTOR_H__
#define __DNAALIGNFEATUREADAPTOR_H__

#include "BaseFeatureAdaptor.h"
#include "AdaptorTypes.h"

struct DNAAlignFeatureAdaptorStruct {
  BASEFEATUREADAPTOR_DATA
};

DNAAlignFeatureAdaptor *DNAAlignFeatureAdaptor_new(DBAdaptor *dba);
int DNAAlignFeatureAdaptor_store(BaseFeatureAdaptor *bfa, Vector *features);
NameTableType *DNAAlignFeatureAdaptor_getTables(void); 
char **DNAAlignFeatureAdaptor_getColumns(void);
NameTableType *DNAAlignFeatureAdaptor_leftJoin(void);
Vector *DNAAlignFeatureAdaptor_objectsFromStatementHandle(BaseFeatureAdaptor *bfa,
                                                       StatementHandle *sth,
                                                       AssemblyMapper *assMapper,
                                                       Slice *slice);


#define DNAAlignFeatureAdaptor_fetchByDbID(dafa, id) BaseFeatureAdaptor_fetchByDbID((BaseFeatureAdaptor *)(dafa), (id))
#define DNAAlignFeatureAdaptor_clearCache(dafa) BaseFeatureAdaptor_clearCache((BaseFeatureAdaptor *)(dafa))
#define DNAAlignFeatureAdaptor_fetchAllBySliceAndScore(dafa, slice, score, lname) \
          BaseFeatureAdaptor_fetchAllBySliceAndScore((BaseFeatureAdaptor *)(dafa), (slice), (score), (lname))
#define DNAAlignFeatureAdaptor_fetchAllByRawContigAndScore(dafa, contig, score, lname) \
          BaseFeatureAdaptor_fetchAllByRawContigAndScore((BaseFeatureAdaptor *)(dafa), (contig), (score), (lname))


#endif
