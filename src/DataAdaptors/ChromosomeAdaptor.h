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

#ifndef __CHROMOSOMEADAPTOR_H__
#define __CHROMOSOMEADAPTOR_H__

#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "Chromosome.h"
#include "IDHash.h"
#include "StringHash.h"

struct ChromosomeAdaptorStruct {
  BASEADAPTOR_DATA
  IDHash *chrCache;
  StringHash *chrNameCache;
};

ChromosomeAdaptor *ChromosomeAdaptor_new(DBAdaptor *dba);
Chromosome *ChromosomeAdaptor_fetchByDbID(ChromosomeAdaptor *ca, IDType dbID);
Chromosome *ChromosomeAdaptor_fetchByChrName(ChromosomeAdaptor *ca, char *chrName);

Chromosome *ChromosomeAdaptor_chromosomeFromRow(ChromosomeAdaptor *ca, ResultRow  *row);


#endif
