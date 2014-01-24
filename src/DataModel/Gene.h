/*
 * Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *      http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef __GENE_H__
#define __GENE_H__

#include "DataModelTypes.h"
#include "AnnotatedSeqFeature.h"
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
  Vector *transcripts;
  ECOSTRING biotype;
  ECOSTRING status;
  ECOSTRING source;
  char *externalName;
  ECOSTRING externalDb;
  ECOSTRING externalStatus;
  char *description;
  char startIsSet;
  char endIsSet;
  char strandIsSet;
  IDType canonicalTranscriptId;
  char *canonicalAnnotation;
  Vector *dbLinks;
  Vector *attributes;
};
#undef FUNCSTRUCTTYPE

Gene *Gene_new(void);

#define Gene_isStored(gene, db) Storable_isStored(&((gene)->st), (db))

#define Gene_setDbID(gene,dbID) AnnotatedSeqFeature_setDbID((gene),dbID)
#define Gene_getDbID(gene) AnnotatedSeqFeature_getDbID((gene))

#define Gene_setStart(gene,start) AnnotatedSeqFeature_setStart((gene),start)
#define Gene_getStart(gene) AnnotatedSeqFeature_getStart((gene))

#define Gene_setEnd(gene,end) AnnotatedSeqFeature_setEnd((gene),end)
#define Gene_getEnd(gene) AnnotatedSeqFeature_getEnd((gene))

#define Gene_setStrand(gene,strand) AnnotatedSeqFeature_setStrand((gene),strand)
#define Gene_getStrand(gene) AnnotatedSeqFeature_getStrand((gene))

#define Gene_getLength(gene) AnnotatedSeqFeature_getLength((gene))

#define Gene_setSlice(gene,sl) AnnotatedSeqFeature_setSlice((gene),(sl))
#define Gene_getSlice(gene) AnnotatedSeqFeature_getSlice((gene))

#define Gene_setAdaptor(gene,ad) AnnotatedSeqFeature_setAdaptor((gene),ad)
#define Gene_getAdaptor(gene) AnnotatedSeqFeature_getAdaptor((gene))

#define Gene_setDisplayXref(gene,xref) AnnotatedSeqFeature_setDisplayXref((gene),xref)
#define Gene_getDisplayXref(gene) AnnotatedSeqFeature_getDisplayXref((gene))

#define Gene_getSeqRegionStart(g) SeqFeature_getSeqRegionStart((g))
#define Gene_getSeqRegionEnd(g) SeqFeature_getSeqRegionEnd((g))
#define Gene_getSeqRegionStrand(g) SeqFeature_getSeqRegionStrand((g))

#define Gene_setStableId(gene,stableId)  StableIdInfo_setStableId(&((gene)->si),stableId)
char *Gene_getStableId(Gene *gene);

ECOSTRING Gene_setBiotype(Gene *gene, char *biotype);
#define Gene_getBiotype(gene)  (gene)->biotype

ECOSTRING Gene_setExternalStatus(Gene *gene, char *externalStatus);
#define Gene_getExternalStatus(gene)  (gene)->externalStatus

ECOSTRING Gene_setExternalDb(Gene *gene, char *externalDb);
#define Gene_getExternalDb(gene)  (gene)->externalDb


ECOSTRING Gene_setStatus(Gene *gene, char *status);
#define Gene_getStatus(gene)  (gene)->status

ECOSTRING Gene_setSource(Gene *gene, char *source);
#define Gene_getSource(gene)  (gene)->source

char *Gene_setDescription(Gene *gene, char *description);
#define Gene_getDescription(gene)  (gene)->description

char *Gene_setExternalName(Gene *gene, char *externalName);
#define Gene_getExternalName(gene)  (gene)->externalName

#define Gene_setCanonicalTranscriptId(gene, id)  (gene)->canonicalTranscriptId = (id)
#define Gene_getCanonicalTranscriptId(gene)  (gene)->canonicalTranscriptId

char *Gene_setCanonicalAnnotation(Gene *g, char *canonicalAnnotation);
#define Gene_getCanonicalAnnotation(gene)  (gene)->canonicalAnnotation

#define Gene_setCreated(gene,cd)  StableIdInfo_setCreated(&((gene)->si),cd)
#define Gene_getCreated(gene)  StableIdInfo_getCreated(&((gene)->si))

#define Gene_setModified(gene,mod)  StableIdInfo_setModified(&((gene)->si),mod)
#define Gene_getModified(gene)  StableIdInfo_getModified(&((gene)->si))

#define Gene_setVersion(gene,ver)  StableIdInfo_setVersion(&((gene)->si),ver)
#define Gene_getVersion(gene)  StableIdInfo_getVersion(&((gene)->si))

#define Gene_setIsCurrent(gene,isC)  StableIdInfo_setIsCurrent(&((gene)->si),(isC))
#define Gene_getIsCurrent(gene)  StableIdInfo_getIsCurrent(&((gene)->si))

#define Gene_setAnalysis(gene,ana) AnnotatedSeqFeature_setAnalysis((gene),ana)
#define Gene_getAnalysis(gene) AnnotatedSeqFeature_getAnalysis((gene))

//#define Gene_addTranscript(gene,trans) Vector_addElement((gene)->transcripts,(trans))
void Gene_addTranscript(Gene *gene, Transcript *trans);
void Gene_recalculateCoordinates(Gene *gene);

Vector *Gene_getAllTranscripts(Gene *gene);

#define Gene_getTranscriptAt(gene,ind) (Transcript *)Vector_getElementAt((gene)->transcripts,ind)

#define Gene_getTranscriptCount(gene) Vector_getNumElement((gene)->transcripts)

#define Gene_EachTranscript(gene,trans,iter) \
    for (iter=0; iter<Gene_getTranscriptCount(gene); iter++) { \
      trans = Gene_getTranscriptAt(gene,iter);

Gene *Gene_transfer(Gene *gene, Slice *slice);
Gene *Gene_transform(Gene *gene, char *csName, char *csVersion, Slice *toSlice);

Vector *Gene_getAllExons(Gene *gene);
int Gene_getExonCount(Gene *gene);

Gene *Gene_transformToSlice(Gene *gene, Slice *slice);
Gene *Gene_transformToRawContig(Gene *gene);

Vector *Gene_getAllAttributes(Gene *gene, char *attribCode);

Vector *Gene_getAllDBLinks(Gene *g);
Vector *Gene_getAllDBEntries(Gene *g);
int Gene_addDBLink(Gene *gene, DBEntry *dbe);

void Gene_free(Gene *gene);
Gene *Gene_shallowCopy(Gene *gene);


#ifdef __GENE_MAIN__
  GeneFuncs 
    geneFuncs = {
                 Gene_free,
                 Gene_shallowCopy, // shallowCopy
                 NULL, // deepCopy
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
