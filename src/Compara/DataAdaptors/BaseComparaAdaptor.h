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

#ifndef __BASECOMPARAADAPTOR_H__
#define __BASECOMPARAADAPTOR_H__

#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "ComparaDBAdaptor.h"

#define BASECOMPARAADAPTOR_DATA BASEADAPTOR_DEF(ComparaDBAdaptor)
struct BaseComparaAdaptorStruct {
  BASECOMPARAADAPTOR_DATA
};

void BaseComparaAdaptor_init(BaseComparaAdaptor *ba, ComparaDBAdaptor *dba, int adaptorType);

#endif
