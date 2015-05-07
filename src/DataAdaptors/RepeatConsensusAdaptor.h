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

#ifndef __REPEATCONSENSUSADAPTOR_H__
#define __REPEATCONSENSUSADAPTOR_H__

#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "RepeatConsensus.h"
#include "Vector.h"

struct RepeatConsensusAdaptorStruct {
  BASEADAPTOR_DATA
};

RepeatConsensusAdaptor *RepeatConsensusAdaptor_new(DBAdaptor *dba);
RepeatConsensus *RepeatConsensusAdaptor_fetchByDbID(RepeatConsensusAdaptor *ca, IDType dbID);
Vector *RepeatConsensusAdaptor_genericFetch(RepeatConsensusAdaptor *rca, char *whereClause);
RepeatConsensus *RepeatConsensusAdaptor_fetchByName(RepeatConsensusAdaptor *rca, char *name);
RepeatConsensus *RepeatConsensusAdaptor_fetchByNameAndClass(RepeatConsensusAdaptor *rca, char *name, char *class);
Vector *RepeatConsensusAdaptor_fetchByClassAndSeq(RepeatConsensusAdaptor *rca, char *class, char *seq);
int RepeatConsensusAdaptor_store(RepeatConsensusAdaptor *rca, Vector *consensi);
int RepeatConsensusAdaptor_free(RepeatConsensusAdaptor *rc);


#endif
