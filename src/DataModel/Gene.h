#ifndef __GENE_H__
#define __GENE_H__

#include "DataModelTypes.h"
#include "AnnotatedSeqFeature.h"
#include "FeatureSet.h"
#include "StableIdInfo.h"
#include "Slice.h"
#include "Transcript.h"
#include "Vector.h"

typedef struct GeneFuncsStruct {
  ANNOTATEDSEQFEATUREFUNCS_DATA
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

#define Gene_setAdaptor(gene,ad) AnnotatedSeqFeature_setAdaptor((gene),ad)
#define Gene_getAdaptor(gene) AnnotatedSeqFeature_getAdaptor((gene))

#define Gene_setStableId(gene,stableId)  StableIdInfo_setStableId(&((gene)->si),stableId)
char *Gene_getStableId(Gene *gene);

char *Gene_setType(Gene *gene, char *type);
#define Gene_getType(gene)  (gene)->type

#define Gene_setStartIsSet(gene, flag)  (gene)->startIsSet = (flag)
#define Gene_getStartIsSet(gene)  (gene)->startIsSet

#define Gene_setEndIsSet(gene, flag)  (gene)->endIsSet = (flag)
#define Gene_getEndIsSet(gene)  (gene)->endIsSet

#define Gene_setStrandIsSet(gene, flag)  (gene)->strandIsSet = (flag)
#define Gene_getStrandIsSet(gene)  (gene)->strandIsSet

#define Gene_setCreated(gene,cd)  StableIdInfo_setCreated(&((gene)->si),cd)
#define Gene_getCreated(gene)  StableIdInfo_getCreated(&((gene)->si))

#define Gene_setModified(gene,mod)  StableIdInfo_setModified(&((gene)->si),mod)
#define Gene_getModified(gene)  StableIdInfo_getModified(&((gene)->si))

#define Gene_setVersion(gene,ver)  StableIdInfo_setVersion(&((gene)->si),ver)
#define Gene_getVersion(gene)  StableIdInfo_getVersion(&((gene)->si))

int Gene_setStart(SeqFeature *gene,int start);
int Gene_getStart(SeqFeature *gene);

#define Gene_setAnalysis(gene,ana) AnnotatedSeqFeature_setAnalysis((gene),ana)
#define Gene_getAnalysis(gene) AnnotatedSeqFeature_getAnalysis((gene))

int Gene_setEnd(SeqFeature *gene,int end);
int Gene_getEnd(SeqFeature *gene);

int Gene_setStrand(SeqFeature *gene,int strand);
int Gene_getStrand(SeqFeature *gene);

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

#ifdef __GENE_MAIN__
  GeneFuncs 
    geneFuncs = {
                 Gene_getStart, 
                 Gene_setStart,
                 Gene_getEnd,
                 Gene_setEnd,
                 Gene_getStrand,
                 Gene_setStrand,
                 NULL, // getSeq
                 NULL, // setSeq
                 NULL, // getLength,
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
