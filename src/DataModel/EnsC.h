#ifndef __ENSC_H__
#define __ENSC_H__

#include <stdio.h>
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

#ifdef __ECOS_MAIN__
 ECOSTRTABLE *ecoSTable;
#else
 extern ECOSTRTABLE *ecoSTable;
#endif

#endif
