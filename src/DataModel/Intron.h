#ifndef __INTRON_H__
#define __INTRON_H__

#include "EnsC.h"
#include "DataModelTypes.h"
#include "SeqFeature.h"
#include "Slice.h"
#include "Exon.h"

SEQFEATUREFUNC_TYPES(Intron)

typedef struct IntronFuncsStruct {
  SEQFEATUREFUNCS_DATA(Intron)
} IntronFuncs;



#define INTRON_DATA \
  SEQFEATURE_DATA \
  Exon *prevExon; \
  Exon *nextExon;

#define FUNCSTRUCTTYPE IntronFuncs
struct IntronStruct {
  INTRON_DATA
};
#undef FUNCSTRUCTTYPE

int Intron_getIsSpliceCanonical(Intron *intron);
long Intron_getLength(Intron *intron);

#define Intron_setPrevExon(intron,exon) (intron)->prevExon = (exon)
#define Intron_getPrevExon(intron) (intron)->prevExon

#define Intron_setNextExon(intron,exon) (intron)->nextExon = (exon)
#define Intron_getNextExon(intron) (intron)->nextExon

#define Intron_setStart(intron,start) SeqFeature_setStart((intron),start)
#define Intron_getStart(intron) SeqFeature_getStart((intron))

#define Intron_setEnd(intron,end) SeqFeature_setEnd((intron),end)
#define Intron_getEnd(intron) SeqFeature_getEnd((intron))

#define Intron_setScore(intron,score) SeqFeature_setScore((intron),score)
#define Intron_getScore(intron) SeqFeature_getScore((intron))

#define Intron_setStrand(intron,strand) SeqFeature_setStrand((intron),(strand))
#define Intron_getStrand(intron) SeqFeature_getStrand((intron))

#define Intron_setAnalysis(intron,ana) SeqFeature_setAnalysis((intron),(ana))
#define Intron_getAnalysis(intron) SeqFeature_getAnalysis((intron))

#define Intron_setSlice(intron,c) SeqFeature_setSlice((intron),(c))
#define Intron_getSlice(intron) SeqFeature_getSlice((intron))

#define Intron_free(intron) SeqFeature_free((intron))

Intron *Intron_new(Exon *e1, Exon *e2, Analysis *analysis);

void Intron_freeImpl(Intron *sf);

#ifdef __INTRON_MAIN__
 IntronFuncs
   intronFuncs = {
                      Intron_freeImpl,
                      NULL, // shallowCopy
                      NULL, // deepCopy
                      NULL, // getStart
                      NULL, // setStart
                      NULL, // getEnd
                      NULL, // setEnd
                      NULL, // getStrand
                      NULL, // setStrand
                      NULL, // getSeq
                      NULL, // setSeq
                      Intron_getLength, // getLength
                      NULL  // reverseComplement
                     };
#else
 extern IntronFuncs intronFuncs;
#endif


#endif
