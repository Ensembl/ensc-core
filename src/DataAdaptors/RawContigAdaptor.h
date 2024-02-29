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

#ifndef __RAWCONTIGADAPTOR_H__
#define __RAWCONTIGADAPTOR_H__

#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "RawContig.h"
#include "IDHash.h"

struct RawContigAdaptorStruct {
  BASEADAPTOR_DATA
  IDHash *rawContigCache;
};

RawContigAdaptor *RawContigAdaptor_new(DBAdaptor *dba);
RawContig *RawContigAdaptor_fetchByDbID(RawContigAdaptor *rca, IDType dbID);
RawContig *RawContigAdaptor_rawContigFromRow(RawContigAdaptor *rca, ResultRow *row);
void RawContigAdaptor_fillRawContigWithRow(RawContigAdaptor *rca, RawContig *rc, ResultRow *row);
void RawContigAdaptor_fetchAttributes(RawContigAdaptor *rca, RawContig *rc);
RawContig *RawContigAdaptor_fetchByName(RawContigAdaptor *rca, char *name);

#endif
