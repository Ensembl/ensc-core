/*
 * Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

#ifndef __SIMPLEFEATUREADAPTOR_H__
#define __SIMPLEFEATUREADAPTOR_H__

#include "BaseFeatureAdaptor.h"
#include "AdaptorTypes.h"
#include "SimpleFeature.h"

struct SimpleFeatureAdaptorStruct {
  BASEFEATUREADAPTOR_DATA
};

SimpleFeatureAdaptor *SimpleFeatureAdaptor_new(DBAdaptor *dba);

int SimpleFeatureAdaptor_store(BaseFeatureAdaptor *baf, Vector *features);
NameTableType *SimpleFeatureAdaptor_getTables(void);
char **SimpleFeatureAdaptor_getColumns(void);
Vector *SimpleFeatureAdaptor_objectsFromStatementHandle(BaseFeatureAdaptor *baf,
                                                     StatementHandle *sth,
                                                     AssemblyMapper *mapper,
                                                     Slice *slice);



#define SimpleFeatureAdaptor_fetchByDbID(sfa, id) BaseFeatureAdaptor_fetchByDbID((BaseFeatureAdaptor *)(sfa), (id))
#define SimpleFeatureAdaptor_fetchAllBySliceAndScore(sfa, slice, score, lname) \
          BaseFeatureAdaptor_fetchAllBySliceAndScore((BaseFeatureAdaptor *)(sfa), (slice), (score), (lname))
#define SimpleFeatureAdaptor_fetchAllByRawContigAndScore(sfa, rc, score, lname) \
          BaseFeatureAdaptor_fetchAllByRawContigAndScore((BaseFeatureAdaptor *)(sfa), (rc), (score), (lname))

#endif
