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

#define __ECOS_MAIN__
#include "EnsC.h"
#undef __ECOS_MAIN__
#include "Stream.h"
#include "StrUtil.h"

void initEnsC(int argc, char **argv) {
  if (!EcoString_initTable(&ecoSTable)) {
    fprintf(stderr, "Failed initialising ecoSTable\n");
    exit(1);
  }

  StrUtil_copyString(&EnsC_progName, argv[0], 0);

  Stream_setDefaults(0);
}


int idTypeCompFunc(const void *one, const void *two) {
  IDType id1 = *((IDType *)one);
  IDType id2 = *((IDType *)two);

  return id1-id2;
}

long *long_new(long val) {
  long *longP;
  if ((longP = (long *)calloc(1,sizeof(long))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for long\n");
    return NULL;
  }
  *longP = val;
  return longP;
}

double *double_new(double val) {
  double *doubleP;
  if ((doubleP = (double *)calloc(1,sizeof(double))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for double\n");
    return NULL;
  }
  *doubleP = val;
  return doubleP;
}

IDType *IDType_new(IDType val) {
  IDType *IDTypeP;
  if ((IDTypeP = (IDType *)calloc(1,sizeof(IDType))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for IDType\n");
    return NULL;
  }
  *IDTypeP = val;
  return IDTypeP;
}

