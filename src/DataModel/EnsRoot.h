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

#ifndef __ENSROOT_H__
#define __ENSROOT_H__

#include "Object.h"

#define ENSROOT_DATA OBJECT_DATA

#define EnsRoot_free(ef) Object_free((ef))

#define EnsRoot_shallowCopy(sf) Object_shallowCopy((sf))
#define EnsRoot_deepCopy(sf) Object_deepCopy((sf))

#endif
