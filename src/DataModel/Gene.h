#ifndef __GENE_H__
#define __GENE_H__

#include "DataModelTypes.h"
#include "AnnotatedSeqFeature.h"
#include "FeatureSet.h"
#include "StableIdInfo.h"
#include "Slice.h"
#include "Transcript.h"
#include "Vector.h"

ANNOTATEDSEQFEATUREFUNC_TYPES(Gene)

typedef struct GeneFuncsStruct {
  ANNOTATEDSEQFEATUREFUNCS_DATA(Gene)
} GeneFuncs;

#define FUNCSTRUCTTYPE GeneFuncs
struct GeneStruct {
  ANNOTATEDSEQFEATURE_DATA
  FeatureSet fs;
  char *type;
  char startIsSet;
  char endIsSet;
  char strandIsSet;
  Vector *dbLinks;
};
#undef FUNCSTRUCTTYPE

Gene *Gene_new(void);

#define Gene_setDbID(gene,dbID) AnnotatedSeqFeature_setDbID((gene),dbID)
#define Gene_getDbID(gene) AnnotatedSeqFeature_getDbID((gene))

#define Gene_setStart(gene,start) AnnotatedSeqFeature_setStart((gene),start)
#define Gene_getStart(gene) AnnotatedSeqFeature_getStart((gene))

#define Gene_setEnd(gene,end) AnnotatedSeqFeature_setEnd((gene),end)
#define Gene_getEnd(gene) AnnotatedSeqFeature_getEnd((gene))

#define Gene_setStrand(gene,strand) AnnotatedSeqFeature_setStrand((gene),strand)
#define Gene_getStrand(gene) AnnotatedSeqFeature_getStrand((gene))

#define Gene_setAdaptor(gene,ad) AnnotatedSeqFeature_setAdaptor((gene),ad)
#define Gene_getAdaptor(gene) AnnotatedSeqFeature_getAdaptor((gene))

#define Gene_setDisplayXref(gene,xref) AnnotatedSeqFeature_setDisplayXref((gene),xref)
#define Gene_getDisplayXref(gene) AnnotatedSeqFeature_getDisplayXref((gene))

#define Gene_setStableId(gene,stableId)  StableIdInfo_setStableId(&((gene)->si),stableId)
char *Gene_getStableId(Gene *gene);

char *Gene_setType(Gene *gene, char *type);
#define Gene_getType(gene)  (gene)->type

#define Gene_setCreated(gene,cd)  StableIdInfo_setCreated(&((gene)->si),cd)
#define Gene_getCreated(gene)  StableIdInfo_getCreated(&((gene)->si))

#define Gene_setModified(gene,mod)  StableIdInfo_setModified(&((gene)->si),mod)
#define Gene_getModified(gene)  StableIdInfo_getModified(&((gene)->si))

#define Gene_setVersion(gene,ver)  StableIdInfo_setVersion(&((gene)->si),ver)
#define Gene_getVersion(gene)  StableIdInfo_getVersion(&((gene)->si))

#define Gene_setAnalysis(gene,ana) AnnotatedSeqFeature_setAnalysis((gene),ana)
#define Gene_getAnalysis(gene) AnnotatedSeqFeature_getAnalysis((gene))

#define Gene_addTranscript(gene,trans) FeatureSet_addFeature(&((gene)->fs),trans)
#define Gene_getTranscriptAt(gene,ind) (Transcript *)FeatureSet_getFeatureAt(&((gene)->fs),ind)

#define Gene_getTranscriptCount(gene) FeatureSet_getNumFeature(&((gene)->fs))

#define Gene_EachTranscript(gene,trans,iter) \
    for (iter=0; iter<Gene_getTranscriptCount(gene); iter++) { \
      trans = Gene_getTranscriptAt(gene,iter);

Gene *Gene_transformToSlice(Gene *gene, Slice *slice);
Vector *Gene_getAllExons(Gene *gene);

Gene *Gene_transformToRawContig(Gene *gene);


Vector *Gene_getAllDBLinks(Gene *g);
int Gene_addDBLink(Gene *gene, DBEntry *dbe);

void Gene_free(Gene *gene);


#ifdef __GENE_MAIN__
  GeneFuncs 
    geneFuncs = {
                 Gene_free,
                 NULL,
                 NULL,
                 NULL,
                 NULL,
                 NULL,
                 NULL,
                 NULL, // getSeq
                 NULL, // setSeq
                 NULL, // getLength,
                 NULL, // reverseComplement
                 Gene_transformToRawContig,
                 Gene_transformToSlice,
                 NULL, // transformRawContigToSlice
                 NULL, // transformSliceToRawContig
                 NULL  // transformSliceToSlice
                };
#else 
  extern GeneFuncs geneFuncs;
#endif

#endif
