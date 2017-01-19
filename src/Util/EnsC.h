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

#ifndef __ENSC_H__
#define __ENSC_H__

#include <stdio.h>
#include <stdlib.h>
#include "EcoString.h"

#ifdef TRUE
#undef TRUE
#endif
#define TRUE 1

#ifdef FALSE
#undef FALSE
#endif
#define FALSE 0


typedef enum CopyDepthEnum {
  NONE_DEPTH,
  SHALLOW_DEPTH,
  DEEP_DEPTH
} CopyDepth;

typedef long long IDType;
#ifdef __osf__
 #define IDFMTSTR "%lld"
#else
 #define IDFMTSTR "%qd"
#endif

#define EXTREMELEN 65536

#define POS_UNDEF   -1111111111L
#define LENGTH_UNDEF POS_UNDEF

#define FLOAT_UNDEF -1111111.0

#define CHARFLAG_UNSET -127
#define PHASE_UNDEF CHARFLAG_UNSET
#define STRAND_UNDEF CHARFLAG_UNSET

#ifdef __ECOS_MAIN__
 ECOSTRTABLE *ecoSTable;
 int trueVal = 1;
 char *EnsC_progName;
 
#else
 extern ECOSTRTABLE *ecoSTable;
 extern int trueVal;
 extern char *EnsC_progName;
#endif

void initEnsC(int argc, char **argv);

typedef int (*SortCompFunc)(const void *a, const void *b);
int idTypeCompFunc(const void *one, const void *two);
long *long_new(long val);
double *double_new(double val);
IDType *IDType_new(IDType val);


#endif
