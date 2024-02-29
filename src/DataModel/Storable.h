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

#ifndef __STORABLE_H__
#define __STORABLE_H__

#include "EnsC.h"
#include "DataModelTypes.h"
#include "AdaptorTypes.h"

struct StorableStruct {
  IDType     dbID;
  BaseAdaptor *adaptor;
};


#define Storable_getDbID(st) (st)->dbID
#define Storable_setDbID(st,id) (st)->dbID = id
#define Storable_getAdaptor(st) (st)->adaptor
#define Storable_setAdaptor(st,ad) (st)->adaptor = ad

int Storable_isStored(Storable *storable, DBAdaptor *db);

#endif
