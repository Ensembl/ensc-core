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

#ifndef __GENOMICALIGNADAPTOR_H__
#define __GENOMICALIGNADAPTOR_H__

#include "BaseComparaAdaptor.h"
#include "ComparaAdaptorTypes.h"
#include "GenomicAlign.h"
#include "Vector.h"

struct GenomicAlignAdaptorStruct {
  BASECOMPARAADAPTOR_DATA
  int maxAlignmentLength;
};

GenomicAlignAdaptor *GenomicAlignAdaptor_new(ComparaDBAdaptor *dba);

Vector *GenomicAlignAdaptor_fetchAllByDNAFragGenomeDB(GenomicAlignAdaptor *gaa,
               DNAFrag *dnaFrag, GenomeDB *targetGenome, int *startP, int *endP,
               char *alignmentType);
Vector *GenomicAlignAdaptor_objectsFromStatementHandle(GenomicAlignAdaptor *gaa, StatementHandle *sth,
                                                       int reverse);
void GenomicAlignAdaptor_addDerivedAlignments(GenomicAlignAdaptor *gaa,
                     Vector *mergedAligns, GenomicAlign *alignA, GenomicAlign *alignB);
IDType GenomicAlignAdaptor_methodLinkIdByAlignmentType(GenomicAlignAdaptor *gaa, char *alignmentType);
char *GenomicAlignAdaptor_alignmentTypeByMethodLinkId(GenomicAlignAdaptor *gaa, IDType methodLinkId);




#endif
