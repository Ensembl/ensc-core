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
