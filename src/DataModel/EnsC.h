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

#define POS_UNDEF -1111111111L
#define LENGTH_UNDEF POS_UNDEF

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


#endif
