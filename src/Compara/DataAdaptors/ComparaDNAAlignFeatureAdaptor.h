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

#ifndef __COMPARADNAALIGNFEATUREADAPTOR_H__
#define __COMPARADNAALIGNFEATUREADAPTOR_H__

#include "BaseComparaAdaptor.h"
#include "ComparaAdaptorTypes.h"
#include "DNAAlignFeature.h"
#include "Vector.h"
#include "Cache.h"

struct ComparaDNAAlignFeatureAdaptorStruct {
  BASECOMPARAADAPTOR_DATA
  Cache *regionCache;
};

ComparaDNAAlignFeatureAdaptor *ComparaDNAAlignFeatureAdaptor_new(ComparaDBAdaptor *dba);
Vector *ComparaDNAAlignFeatureAdaptor_fetchAllBySpeciesRegion(ComparaDNAAlignFeatureAdaptor *dafa,
                                                   char *csSpecies, char *csAssembly,
                                                   char *qySpecies, char *qyAssembly,
                                                   char *chrName, int start, int end,
                                                   char *alignmentType);
Vector *ComparaDNAAlignFeatureAdaptor_fetchAllBySlice(ComparaDNAAlignFeatureAdaptor *dafa,
              Slice *slice, char *qySpecies, char *qyAssembly, char *assemblyType);

#endif
