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

#ifndef __HOMOLOGYADAPTOR_H__
#define __HOMOLOGYADAPTOR_H__

#include "BaseComparaAdaptor.h"
#include "ComparaAdaptorTypes.h"
#include "Homology.h"
#include "Vector.h"

struct HomologyAdaptorStruct {
  BASECOMPARAADAPTOR_DATA
};

HomologyAdaptor *HomologyAdaptor_new(ComparaDBAdaptor *dba);

Vector *HomologyAdaptor_fetchHomologuesOfGeneInSpecies(HomologyAdaptor *ha,
                          char *sp, char *gene, char *hSp);
Vector *HomologyAdaptor_fetchHomologuesBySpeciesRelationshipId(HomologyAdaptor *ha,
                   char *hSpecies, IDType internalId);

IDType HomologyAdaptor_getRelationship(HomologyAdaptor *ha, char *qStr);
int HomologyAdaptor_getRelationships(HomologyAdaptor *ha, char *qStr, IDType **idsP);
Vector *HomologyAdaptor_getHomologues(HomologyAdaptor *ha, char *qStr);
Vector *HomologyAdaptor_listStableIdsFromSpecies(HomologyAdaptor *ha, char *sp);



#endif
