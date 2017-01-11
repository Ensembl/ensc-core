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

#ifndef __RESULTROW_H__
#define __RESULTROW_H__

#include "EnsC.h"
#include "Object.h"

typedef struct ResultRowStruct ResultRow;

typedef char *    (*ResultRow_getStringAtFunc)(ResultRow *row, int ind);
typedef char *    (*ResultRow_getStringAllowNullAtFunc)(ResultRow *row, int ind);
typedef char *    (*ResultRow_getStringCopyAtFunc)(ResultRow *row, int ind);
typedef int       (*ResultRow_getIntAtFunc)(ResultRow *row, int ind);
typedef long      (*ResultRow_getLongAtFunc)(ResultRow *row, int ind);
typedef IDType     (*ResultRow_getLongLongAtFunc)(ResultRow *row, int ind);
typedef double    (*ResultRow_getDoubleAtFunc)(ResultRow *row, int ind);
typedef char *    (*ResultRow_colFunc)(ResultRow *row,int ind);

OBJECTFUNC_TYPES(ResultRow)

typedef struct ResultRowFuncsStruct {
  OBJECTFUNCS_DATA(ResultRow)
} ResultRowFuncs;


#define RESULTROW_DATA \
  OBJECT_DATA \
  ResultRow_getStringAtFunc   getStringAt; \
  ResultRow_getStringAtFunc   getStringAllowNullAt; \
  ResultRow_getStringCopyAtFunc   getStringCopyAt; \
  ResultRow_getIntAtFunc      getIntAt; \
  ResultRow_getLongAtFunc     getLongAt; \
  ResultRow_getLongLongAtFunc getLongLongAt; \
  ResultRow_getDoubleAtFunc   getDoubleAt; \
  ResultRow_colFunc           col;

#define FUNCSTRUCTTYPE ResultRowFuncs
struct ResultRowStruct {
  RESULTROW_DATA
};
#undef FUNCSTRUCTTYPE


#endif
