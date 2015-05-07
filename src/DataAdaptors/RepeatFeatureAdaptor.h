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

#ifndef __REPEATFEATUREADAPTOR_H__
#define __REPEATFEATUREADAPTOR_H__

#include "BaseFeatureAdaptor.h"
#include "AdaptorTypes.h"

struct RepeatFeatureAdaptorStruct {
  BASEFEATUREADAPTOR_DATA
};

RepeatFeatureAdaptor *RepeatFeatureAdaptor_new(DBAdaptor *dba);

int RepeatFeatureAdaptor_store(BaseFeatureAdaptor *bfa, Vector *features);
char *RepeatFeatureAdaptor_defaultWhereClause();
NameTableType *RepeatFeatureAdaptor_getTables();
char **RepeatFeatureAdaptor_getColumns();
Vector *RepeatFeatureAdaptor_fetchAllBySlice(RepeatFeatureAdaptor *rfa, Slice *slice, char *logicName, Vector *repeatTypes);
Vector *RepeatFeatureAdaptor_objectsFromStatementHandle(BaseFeatureAdaptor *bfa,
                                                     StatementHandle *sth,
                                                     AssemblyMapper *mapper,
                                                     Slice *slice);


#define RepeatFeatureAdaptor_fetchByDbID(rfa, id) BaseFeatureAdaptor_fetchByDbID((BaseFeatureAdaptor *)(rfa), (id))
#define RepeatFeatureAdaptor_fetchAllBySliceConstraint(rfa, slice, constraint, lname) \
          BaseFeatureAdaptor_fetchAllBySliceConstraint((BaseFeatureAdaptor *)(rfa), (slice), (constraint), (lname))
#define RepeatFeatureAdaptor_fetchAllByRawContig(rfa, contig, lname) \
          BaseFeatureAdaptor_fetchAllByRawContig((BaseFeatureAdaptor *)(rfa), (contig), (lname))

#endif
