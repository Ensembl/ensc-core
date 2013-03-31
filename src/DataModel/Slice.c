#define __SLICE_MAIN__
#include "Slice.h"
#undef __SLICE_MAIN__

#include "ChromosomeAdaptor.h"
#include "DNAAlignFeatureAdaptor.h"
#include "GeneAdaptor.h"
#include "PredictionTranscriptAdaptor.h"
#include "ProteinAlignFeatureAdaptor.h"
#include "RepeatFeatureAdaptor.h"
#include "SimpleFeatureAdaptor.h"
#include "SliceAdaptor.h"

#include "Gene.h"

#include "Vector.h"
#include "StrUtil.h"

Slice *Slice_new(char *chr, int start, int end, int strand, char *assemblyType,
                 SliceAdaptor *sa, IDType dbID, int empty) {
  Slice *slice;

  if ((slice = (Slice *)calloc(1,sizeof(Slice))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for slice\n");
    return NULL;
  }

  slice->objectType = CLASS_SLICE;
  Object_incRefCount(slice);

  slice->funcs = &sliceFuncs;

  if (!empty) {
    if( !chr || !assemblyType) {
      fprintf(stderr,"ERROR: Do not have all the parameters for slice\n");
      exit(1);
    }
    Slice_setChrName(slice,chr);
    Slice_setChrStart(slice,start);
    Slice_setChrEnd(slice,end);
    Slice_setStrand(slice,strand);
  } else {
    Slice_setStrand(slice,1);
    Slice_setChrStart(slice,1);
    Slice_setEmptyFlag(slice,TRUE);

    // empty Slices are used to do mapping to chromosomal coords.
    // After the mapping, the Slice contains chr_name and is reference
    // point for the mapped object
  }

  Slice_setAssemblyType(slice,assemblyType);
  Slice_setAdaptor(slice,(BaseAdaptor *)sa);
  Slice_setDbID(slice,dbID);

/*
  if (sa && chr) {
    ChromosomeAdaptor *ca = DBAdaptor_getChromosomeAdaptor(sa->dba);
    Chromosome *chromosome = ChromosomeAdaptor_fetchByChrName(ca, chr);
    Slice_setChrId(slice,Chromosome_getDbID(chromosome));
  }
*/

/*
  if( defined $adaptor && !defined $type ) {
    $self->assembly_type
      ( $adaptor->db()->get_MetaContainer()->get_default_assembly());
  }
*/
  return slice;
}

ECOSTRING Slice_getName(Slice *slice) {
  char tmpStr[MAXSTRLEN];

  // NIY printf("Need to redo this\n");
  if (!slice->name) {
    sprintf(tmpStr,"%s.%d-%d",Slice_getChrName(slice),
            Slice_getChrStart(slice),Slice_getChrEnd(slice));
    EcoString_copyStr(ecoSTable,&(slice->name),tmpStr,0);
  }
  return slice->name;
}

/* logicName is optional - set to NULL if you want them all */
Vector *Slice_getAllGenes(Slice *slice, char *logicName) {
  GeneAdaptor *ga = DBAdaptor_getGeneAdaptor(Slice_getAdaptor(slice)->dba);

  return GeneAdaptor_fetchAllBySlice(ga,slice,logicName);
}

Vector *Slice_getAllGenesByType(Slice *slice, char *type) {
  Vector *allGenes = Slice_getAllGenes(slice, NULL);
  Vector *typedGenes = Vector_new();
  int i;

  for (i=0;i<Vector_getNumElement(allGenes);i++) {
    Gene *gene = Vector_getElementAt(allGenes,i);
    if (!strcmp(Gene_getType(gene),type)) {
      Vector_addElement(typedGenes,gene);
    }
  }
  //Vector_free(allGenes,NULL);
  
  return typedGenes;
}

Vector *Slice_getAllSimpleFeatures(Slice *slice, char *logicName, double *score) {

  SimpleFeatureAdaptor *sfa = DBAdaptor_getSimpleFeatureAdaptor(Slice_getAdaptor(slice)->dba);
  return SimpleFeatureAdaptor_fetchAllBySliceAndScore(sfa, slice, score, logicName);
}

Vector *Slice_getAllDNAAlignFeatures(Slice *slice, char *logicName, double *score) {

  DNAAlignFeatureAdaptor *dafa = DBAdaptor_getDNAAlignFeatureAdaptor(Slice_getAdaptor(slice)->dba);
  return DNAAlignFeatureAdaptor_fetchAllBySliceAndScore(dafa, slice, score, logicName);
}

Vector *Slice_getAllDNAPepAlignFeatures(Slice *slice, char *logicName, double *score) {

  ProteinAlignFeatureAdaptor *pafa = DBAdaptor_getProteinAlignFeatureAdaptor(Slice_getAdaptor(slice)->dba);
  return ProteinAlignFeatureAdaptor_fetchAllBySliceAndScore(pafa, slice, score, logicName);
}

Vector *Slice_getAllRepeatFeatures(Slice *slice, char *logicName) {

  RepeatFeatureAdaptor *rfa = DBAdaptor_getRepeatFeatureAdaptor(Slice_getAdaptor(slice)->dba);
  return RepeatFeatureAdaptor_fetchAllBySlice(rfa, slice, logicName);
}

Vector *Slice_getAllPredictionTranscripts(Slice *slice, char *logicName) {

  PredictionTranscriptAdaptor *pta = DBAdaptor_getPredictionTranscriptAdaptor(Slice_getAdaptor(slice)->dba);
  return PredictionTranscriptAdaptor_fetchAllBySlice(pta, slice);
}


ECOSTRING Slice_setChrName(Slice *sl, char *chrName) {
  EcoString_copyStr(ecoSTable, &(sl->chrName), chrName, 0);

  return sl->chrName;
}

ECOSTRING Slice_setAssemblyType(Slice *sl, char *assemblyType) {
  EcoString_copyStr(ecoSTable, &(sl->assemblyType), assemblyType, 0);

  return sl->assemblyType;
}

char *Slice_getSeq(Slice *slice) {
  SequenceAdaptor *sa;

  sa = DBAdaptor_getSequenceAdaptor(Slice_getAdaptor(slice)->dba);

  return SequenceAdaptor_fetchBySliceStartEndStrand(sa,slice,1,-1,1);
}


char *Slice_getSubSeq(Slice *slice, int start, int end, int strand) {
  SequenceAdaptor *sa;

  if (end < start) {
    fprintf(stderr,"Error: End coord is less then start coord\n");
    exit(1);
  }

  if (strand != -1 && strand != 1 ) {
    fprintf(stderr,"Error: Invalid strand [%d] in call to Slice subseq.\n",strand);
    exit(1);
  }

  sa = DBAdaptor_getSequenceAdaptor(Slice_getAdaptor(slice)->dba);
  return SequenceAdaptor_fetchBySliceStartEndStrand(sa,slice,start,end,strand);
}

void Slice_free(Slice *slice) {
  Object_decRefCount(slice);

  if (Object_getRefCount(slice) > 0) {
    return;
  } else if (Object_getRefCount(slice) < 0) {
    fprintf(stderr,"Error: Negative reference count for Slice\n"
                   "       Freeing it anyway\n");
  }

  BaseContig_freePtrs(slice);

  if (slice->name) EcoString_freeStr(ecoSTable, slice->name);
  if (slice->chrName) EcoString_freeStr(ecoSTable, slice->chrName);
  if (slice->assemblyType) EcoString_freeStr(ecoSTable, slice->assemblyType);


  free(slice);
}
