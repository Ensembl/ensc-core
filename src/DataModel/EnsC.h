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
#define STRAND_UNDEF -1231231231L

#define CHARFLAG_UNSET -127

#ifdef __ECOS_MAIN__
 ECOSTRTABLE *ecoSTable;
 int trueVal = 1;
#else
 extern ECOSTRTABLE *ecoSTable;
 extern int trueVal;
#endif

void initEnsC(void);

typedef int (*SortCompFunc)(const void *a, const void *b);
int idTypeCompFunc(const void *one, const void *two);


#endif
