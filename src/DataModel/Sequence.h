#ifndef __SEQUENCE_H__
#define __SEQUENCE_H__

#include "DataModelTypes.h"
#include "EnsRoot.h"


#define SEQUENCEFUNCS_DATA \
  int i; // No funcs

typedef struct SequenceFuncsStruct {
  SEQUENCEFUNCS_DATA
} SequenceFuncs;

#define SEQUENCE_DATA \
  ENSROOT_DATA \
  char *seq; \
  char *name; \
  int length;

#define FUNCSTRUCTTYPE SequenceFuncs
struct SequenceStruct {
  SEQUENCE_DATA
};
#undef FUNCSTRUCTTYPE


#ifdef __SEQUENCE_MAIN__
 SequenceFuncs sequenceFuncs;
#else
 extern SequenceFuncs sequenceFuncs;
#endif


#endif
