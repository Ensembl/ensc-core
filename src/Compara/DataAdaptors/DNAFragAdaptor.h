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

#ifndef __DNAFRAGADAPTOR_H__
#define __DNAFRAGADAPTOR_H__

#include "BaseComparaAdaptor.h"
#include "ComparaAdaptorTypes.h"
#include "DNAFrag.h"
#include "Vector.h"

struct DNAFragAdaptorStruct {
  BASECOMPARAADAPTOR_DATA
};

DNAFragAdaptor *DNAFragAdaptor_new(ComparaDBAdaptor *dba);
DNAFrag *DNAFragAdaptor_fetchByDbID(DNAFragAdaptor *dfa, IDType dbID);
Vector *DNAFragAdaptor_fetchAllByGenomeDBRegion(DNAFragAdaptor *dfa, GenomeDB *genomeDB,
                      char *dnaFragType, char *name, int *startP, int *endP);
Vector *DNAFragAdaptor_fetchAll(DNAFragAdaptor *dfa);
Vector *DNAFragAdaptor_objectsFromStatementHandle(DNAFragAdaptor *dfa, StatementHandle *sth);
IDType DNAFragAdaptor_store(DNAFragAdaptor *dfa, DNAFrag *dnaFrag);
IDType DNAFragAdaptor_storeIfNeeded(DNAFragAdaptor *dfa, DNAFrag *dnaFrag);


#endif
