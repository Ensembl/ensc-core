#ifndef __SEQUENCE_H__
#define __SEQUENCE_H__

#include "DataModelTypes.h"
#include "EnsRoot.h"
#include "EnsC.h"

#define SEQUENCEFUNC_TYPES(CLASSTYPE) \
  typedef void (*CLASSTYPE ## _FreeFunc)(CLASSTYPE *seq); \
  typedef ECOSTRING (*CLASSTYPE ## _GetNameFunc)(CLASSTYPE *seq); \
  typedef char * (*CLASSTYPE ## _GetSeqFunc)(CLASSTYPE *seq); \
  typedef char * (*CLASSTYPE ## _GetSubSeqFunc)(CLASSTYPE *seq, int start, int end, int strand);


#define SEQUENCEFUNCS_DATA(CLASSTYPE) \
  CLASSTYPE ## _FreeFunc free; \
  CLASSTYPE ## _GetNameFunc getName; \
  CLASSTYPE ## _GetSeqFunc getSeq; \
  CLASSTYPE ## _GetSubSeqFunc getSubSeq;

SEQUENCEFUNC_TYPES(Sequence)

typedef struct SequenceFuncsStruct {
  SEQUENCEFUNCS_DATA(Sequence)
} SequenceFuncs;

#define SEQUENCE_DATA \
  ENSROOT_DATA \
  char *seq; \
  ECOSTRING name; \
  int length;

#define FUNCSTRUCTTYPE SequenceFuncs
struct SequenceStruct {
  SEQUENCE_DATA
};
#undef FUNCSTRUCTTYPE

void Sequence_freePtrs(Sequence *seq);
Sequence *Sequence_new();

#ifdef __SEQUENCE_MAIN__
 SequenceFuncs sequenceFuncs;
#else
 extern SequenceFuncs sequenceFuncs;
#endif


#endif
