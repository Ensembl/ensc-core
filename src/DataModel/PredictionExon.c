#define __PREDICTIONEXON_MAIN__
#include "PredictionExon.h"
#undef __PREDICTIONEXON__MAIN__
#include "StrUtil.h"

#include <string.h>

PredictionExon *PredictionExon_new() {
  PredictionExon *pe;

  if ((pe = (PredictionExon *)calloc(1,sizeof(PredictionExon))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for pe\n");
    return NULL;
  }

  pe->objectType = CLASS_PREDICTIONEXON;
  Object_incRefCount(pe);

  pe->funcs = &predictionExonFuncs;
 
  return pe;
}

ECOSTRING PredictionExon_setDisplayLabel(PredictionExon *pe, char *label) {
  EcoString_copyStr(ecoSTable, &(pe->displayLabel), label, 0);
  return pe->displayLabel;
}

PredictionExon *PredictionExon_shallowCopyImpl(PredictionExon *exon) {
  PredictionExon *newPredictionExon = PredictionExon_new();

  memcpy(newPredictionExon,exon,sizeof(PredictionExon));

  return newPredictionExon;
}


void PredictionExon_freeImpl(PredictionExon *pe) {
  Object_decRefCount(pe);

  if (Object_getRefCount(pe) > 0) {
    return;
  } else if (Object_getRefCount(pe) < 0) {
    fprintf(stderr,"Error: Negative reference count for PredictionExon\n"
                   "       Freeing it anyway\n");
  }

  if (pe->displayLabel) EcoString_freeStr(ecoSTable, pe->displayLabel);

  SeqFeature_freePtrs((SeqFeature *)pe);
  free(pe);
}

char  *PredictionExon_getSeqString(PredictionExon *exon) {
  char *seq;

  if (PredictionExon_getSeqCacheString(exon)) {
    return PredictionExon_getSeqCacheString(exon);
  }

  if (!PredictionExon_getSlice(exon)) {
    fprintf(stderr, "Warning: this exon %s doesn't have a contig you won't get a seq\n", PredictionExon_getDisplayLabel(exon));
    return NULL;
  } else {

    seq = BaseContig_getSubSeq(PredictionExon_getSlice(exon),
                               PredictionExon_getStart(exon),
                               PredictionExon_getEnd(exon),
                               1);

    if (PredictionExon_getStrand(exon) == -1){
      SeqUtil_reverseComplement(seq,strlen(seq));
    }
  }
  PredictionExon_setSeqCacheString(exon, seq);

  return PredictionExon_getSeqCacheString(exon);
}

void PredictionExon_loadGenomicMapper(Exon *exon, Mapper *mapper, IDType id, int start) {

// NIY Make the Exon_getContig consistent
  Mapper_addMapCoordinates( mapper, id, start, start+Exon_getLength(exon)-1,
                            Exon_getStrand(exon), Slice_getSeqRegionId(PredictionExon_getSlice(exon)),
                            Exon_getStart(exon),  Exon_getEnd(exon) );
}
