#ifndef __GENE_H__
#define __GENE_H__

#include "DataModelTypes.h"
#include "SeqFeature.h"
#include "FeatureSet.h"
#include "StableIdInfo.h"
#include "Slice.h"
#include "Transcript.h"
#include "Set.h"

struct GeneStruct {
  FeatureSet fs;
  SeqFeature sf;
  StableIdInfo si;
  char *type;
};

Gene *Gene_new(void);

#define Gene_setDbID(g,dbID) SeqFeature_setDbID(&((g)->sf),dbID)
#define Gene_getDbID(g) SeqFeature_getDbID(&((g)->sf))

#define Gene_setAdaptor(g,ad) SeqFeature_setAdaptor(&((g)->sf),ad)
#define Gene_getAdaptor(g) SeqFeature_getAdaptor(&((g)->sf))

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

#define Gene_setStart(gene,start) SeqFeature_setStart(&((gene)->sf),start)
#define Gene_getStart(gene) SeqFeature_getStart(&((gene)->sf))

#define Gene_setAnalysis(gene,ana) SeqFeature_setAnalysis(&((gene)->sf),ana)
#define Gene_getAnalysis(gene) SeqFeature_getAnalysis(&((gene)->sf))

#define Gene_setEnd(gene,end) SeqFeature_setEnd(&((gene)->sf),end)
#define Gene_getEnd(gene) SeqFeature_getEnd(&((gene)->sf))

#define Gene_setStrand(gene,strand) SeqFeature_setStrand(&((gene)->sf),strand)
#define Gene_getStrand(gene) SeqFeature_getStrand(&((gene)->sf))

#define Gene_addTranscript(gene,trans) FeatureSet_addFeature(&((gene)->fs),trans)
#define Gene_getTranscriptAt(gene,ind) (Transcript *)FeatureSet_getFeatureAt(&((gene)->fs),ind)

#define Gene_getTranscriptCount(gene) FeatureSet_getNumFeature(&((gene)->fs))

#define Gene_EachTranscript(gene,trans,iter) \
    for (iter=0; iter<Gene_getTranscriptCount(gene); iter++) { \
      trans = Gene_getTranscriptAt(gene,iter);

Gene *Gene_transformToSlice(Gene *gene, Slice *slice);
Set *Gene_getAllExons(Gene *gene);



#endif
