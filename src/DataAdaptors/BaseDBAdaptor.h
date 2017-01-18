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

#ifndef __BASEDBADAPTOR_H__
#define __BASEDBADAPTOR_H__

#include "DBConnection.h"
#include "AdaptorTypes.h"
#include "EnsC.h"

#define BASEDBADAPTOR_DATA \
  DBConnection  *dbc; \
  MetaContainer *metaContainer; \
  MetaCoordContainer *metaCoordContainer;

struct BaseDBAdaptorStruct {
  BASEDBADAPTOR_DATA
};


#define BaseDBAdaptor_prepare(dba,qStr,qLen) (dba)->dbc->prepare((dba)->dbc,(qStr),(qLen))


#endif
