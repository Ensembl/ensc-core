#ifndef __RESULTROW_H__
#define __RESULTROW_H__

#include "EnsC.h"
#include "Object.h"

typedef struct ResultRowStruct ResultRow;

typedef char *    (*ResultRow_getStringAtFunc)(ResultRow *row, int ind);
typedef int       (*ResultRow_getIntAtFunc)(ResultRow *row, int ind);
typedef long      (*ResultRow_getLongAtFunc)(ResultRow *row, int ind);
typedef int64     (*ResultRow_getLongLongAtFunc)(ResultRow *row, int ind);
typedef double    (*ResultRow_getDoubleAtFunc)(ResultRow *row, int ind);

#define RESULTROW_DATA \
  OBJECT_DATA \
  ResultRow_getStringAtFunc   getStringAt; \
  ResultRow_getIntAtFunc      getIntAt; \
  ResultRow_getLongAtFunc     getLongAt; \
  ResultRow_getLongLongAtFunc getLongLongAt; \
  ResultRow_getDoubleAtFunc   getDoubleAt;

struct ResultRowStruct {
  RESULTROW_DATA
};

#endif
