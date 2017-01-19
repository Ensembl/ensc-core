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

#ifndef __GENOMEDBADAPTOR_H__
#define __GENOMEDBADAPTOR_H__

#include "BaseComparaAdaptor.h"
#include "ComparaAdaptorTypes.h"
#include "IDHash.h"
#include "StringHash.h"
#include "GenomeDB.h"
#include "Vector.h"

struct GenomeDBAdaptorStruct {
  BASECOMPARAADAPTOR_DATA
  IDHash *genomeDBCache;
  StringHash *genomeConsensusXrefList;
  StringHash *genomeQueryXrefList;
};

GenomeDBAdaptor *GenomeDBAdaptor_new(ComparaDBAdaptor *dba);
GenomeDB *GenomeDBAdaptor_fetchByDbID(GenomeDBAdaptor *gda, IDType dbID);
Vector *GenomeDBAdaptor_fetchAll(GenomeDBAdaptor *gda);
GenomeDB *GenomeDBAdaptor_fetchByNameAssembly(GenomeDBAdaptor *gda, char *name, char *assembly);
IDType GenomeDBAdaptor_store(GenomeDBAdaptor *gda, GenomeDB *gdb);
void GenomeDBAdaptor_createGenomeDBs(GenomeDBAdaptor *gda);
int GenomeDBAdaptor_checkForConsensusDb(GenomeDBAdaptor *gda, GenomeDB *queryGdb,
                                        GenomeDB *conGdb, IDType methodLinkId);
int GenomeDBAdaptor_checkForQueryDb(GenomeDBAdaptor *gda, GenomeDB *conGdb,
                                    GenomeDB *queryGdb, IDType methodLinkId);
Vector *GenomeDBAdaptor_getAllDbLinks(GenomeDBAdaptor *gda, GenomeDB *refGdb, IDType methodLinkId);







#endif
