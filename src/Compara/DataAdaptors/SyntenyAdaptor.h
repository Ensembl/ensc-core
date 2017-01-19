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

#ifndef __SYNTEYADAPTOR_H__
#define __SYNTEYADAPTOR_H__

#include "BaseComparaAdaptor.h"
#include "ComparaAdaptorTypes.h"
#include "Vector.h"
#include "SyntenyRegion.h"

struct SyntenyAdaptorStruct {
  BASECOMPARAADAPTOR_DATA
  char *speciesMain;
  char *speciesSecondary;
};

SyntenyAdaptor *SyntenyAdaptor_new(ComparaDBAdaptor *dba);

#define SyntenyAdaptor_setSpeciesMain(sa,sm) (sa)->speciesMain = (sm)
#define SyntenyAdaptor_getSpeciesMain(sa) (sa)->speciesMain

#define SyntenyAdaptor_setSpeciesSecondary(sa,ss) (sa)->speciesSecondary = (ss)
#define SyntenyAdaptor_getSpeciesSecondary(sa) (sa)->speciesSecondary

void SyntenyAdaptor_setSpecies(SyntenyAdaptor *sa, char *species1, char *species2);
Vector *SyntenyAdaptor_getSyntenyForChromosome(SyntenyAdaptor *sa, char *chr, int *startP, int *endP);

#endif
