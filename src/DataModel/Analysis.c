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

#define __ANALYSIS_MAIN__
#include "Analysis.h"
#undef __ANALYSIS_MAIN__

#include <string.h>

Analysis *Analysis_new() {
  Analysis *anal;

  if ((anal = (Analysis *)calloc(1,sizeof(Analysis))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for anal\n");
    return NULL;
  }

  anal->objectType = CLASS_ANALYSIS;

  anal->funcs = &analysisFuncs;

  Object_incRefCount(anal);

  return anal;
}

#define COMPARE_STRINGS(func,a,b) \
  str1 = func((a)); \
  str2 = func((b)); \
  \
  if (str1 && !str2) { \
    return 1; \
  } else if (!str1 && str2) { \
    return -1; \
  } else if (str1 && str2) { \
    if (strcmp(str1,str2)) { return -1; } \
  }

#define COMPARE_INTS(func,a,b) \
  i1 = func((a)); \
  i2 = func((b)); \
  \
  if (i1 && !i2) { \
    return 1; \
  } else if (!i1 && i2) { \
    return -1; \
  } else if (i1 && i2) { \
    if (i1 != i2) { return -1; } \
  }

int Analysis_compare(Analysis *a, Analysis *b) {
  char *str1,*str2;
  int i1,i2;

  COMPARE_STRINGS(Analysis_getLogicName,a,b);
  COMPARE_STRINGS(Analysis_getProgram,a,b);
  COMPARE_INTS(Analysis_getProgramVersion,a,b);
  COMPARE_STRINGS(Analysis_getProgramFile,a,b);
  COMPARE_STRINGS(Analysis_getDb,a,b);
  COMPARE_INTS(Analysis_getDbVersion,a,b);
  COMPARE_STRINGS(Analysis_getDbFile,a,b);
  COMPARE_STRINGS(Analysis_getGFFSource,a,b);
  COMPARE_STRINGS(Analysis_getGFFFeature,a,b);
  COMPARE_STRINGS(Analysis_getModule,a,b);
  COMPARE_INTS(Analysis_getModuleVersion,a,b);
  COMPARE_STRINGS(Analysis_getParameters,a,b);
  
  return 0;
}

void Analysis_free(Analysis *anal) {
  Object_decRefCount(anal);

  if (Object_getRefCount(anal) > 0) {
    return;
  } else if (Object_getRefCount(anal) < 0) {
    fprintf(stderr,"Error: Negative reference count for Analysis\n"
                   "       Freeing it anyway\n");
  }

  if (anal->db)          free(anal->db);
  if (anal->dbFile)      free(anal->dbFile);
  if (anal->program)     free(anal->program);
  if (anal->programFile) free(anal->programFile);
  if (anal->gffSource)   free(anal->gffSource);
  if (anal->gffFeature)  free(anal->gffFeature);
  if (anal->module)      free(anal->module);
  if (anal->parameters)  free(anal->parameters);
  if (anal->created)     free(anal->created);
  if (anal->logicName)   free(anal->logicName);

  free(anal);
}

