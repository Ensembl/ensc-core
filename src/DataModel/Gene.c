#define __GENE_MAIN__
#include "Gene.h"
#undef __GENE_MAIN__

#include "DBAdaptor.h"
#include "GeneAdaptor.h"
#include "IDHash.h"

#include "DBEntryAdaptor.h"

Gene *Gene_new() {
  Gene *gene;

  if ((gene = (Gene *)calloc(1,sizeof(Gene))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for gene\n");
    return NULL;
  }

  Gene_setModified(gene,0);
  Gene_setCreated(gene,0);
  Gene_setVersion(gene,-1);

  gene->objectType = CLASS_GENE;
  Object_incRefCount(gene);

  gene->funcs = &geneFuncs;

  return gene;
}

Vector *Gene_getAllDBLinks(Gene *g) {
  if (!g->dbLinks) {
    GeneAdaptor *ga = (GeneAdaptor *)Gene_getAdaptor(g);

    if (ga) {
      DBEntryAdaptor *dbea = DBAdaptor_getDBEntryAdaptor(ga->dba);
      DBEntryAdaptor_fetchAllByGene(dbea,g);
    } else {
      g->dbLinks = emptyVector;
    }
  }

  return g->dbLinks;
}

int Gene_addDBLink(Gene *g, DBEntry *dbe) {
  if (!g->dbLinks) {
    g->dbLinks = Vector_new();
  }

  Vector_addElement(g->dbLinks, dbe); 
  return 1;
}

char *Gene_getStableId(Gene *gene) {
  GeneAdaptor *ga = (GeneAdaptor *)Gene_getAdaptor(gene);

  if (StableIdInfo_getStableId(&(gene->si)) == NULL && ga) {
    GeneAdaptor_getStableEntryInfo(ga,gene);
  }
  return StableIdInfo_getStableId(&(gene->si));
}

char *Gene_setType(Gene *g, char *type) {
  if ((g->type = (char *)malloc(strlen(type)+1)) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for gene type\n");
    return NULL;
  }

  strcpy(g->type,type);

  return g->type;
}


Vector *Gene_getAllExons(Gene *gene) {
  IDHash *exonHash = IDHash_new(IDHASH_SMALL);
  int i;
  Vector *exonVector = Vector_new();
  void **values;

  for (i=0;i<Gene_getTranscriptCount(gene);i++) {
    Transcript *trans = Gene_getTranscriptAt(gene,i);
    int j;

    for (j=0;j<Transcript_getExonCount(trans);j++) {
      Exon *exon = Transcript_getExonAt(trans,j);
      if (!IDHash_contains(exonHash,(long)exon)) {
        IDHash_add(exonHash,(long)exon,exon);
      }
    }
  }

  values = IDHash_getValues(exonHash);

  for (i=0;i<IDHash_getNumValues(exonHash);i++) {
    Vector_addElement(exonVector, values[i]);
  }
  free(values);
  IDHash_free(exonHash,NULL);
  
  return exonVector;
}

Gene *Gene_transformToSlice(Gene *gene, Slice *slice) {
  IDHash *exonTransforms = IDHash_new(IDHASH_SMALL);
  int i;
  Vector *exons = Gene_getAllExons(gene);

  // transform Exons
  for (i=0;i<Vector_getNumElement(exons); i++) {
    Exon *exon = (Exon *)Vector_getElementAt(exons,i);
     
    Exon *newExon = Exon_transformToSlice(exon,slice);
    IDHash_add(exonTransforms, (IDType)exon, newExon);
  }

  // now need to re-jiggle the transcripts and their
  // translations to account for the re-mapping process

  for (i=0;i<Gene_getTranscriptCount(gene);i++) {
    Transcript *trans = Gene_getTranscriptAt(gene,i);

    // need to grab the translation before starting to
    // re-jiggle the exons

    Transcript_transform(trans, exonTransforms);
  }

  // unset the start, end, and strand - they need to be recalculated

#ifdef DONE
  $self->{_chr_name} = undef;
#endif

  IDHash_free(exonTransforms, NULL);

  return gene;
}

Gene *Gene_transformToRawContig(Gene *gene) {
  IDHash *exonTransforms = IDHash_new(IDHASH_SMALL);
  int i;
  Vector *exons = Gene_getAllExons(gene);

  // transform Exons
  for (i=0;i<Vector_getNumElement(exons); i++) {
    Exon *exon = Vector_getElementAt(exons,i);
     
    Exon *newExon = Exon_transformToRawContig(exon);
    IDHash_add(exonTransforms, (IDType)exon, newExon);
  }

  // now need to re-jiggle the transcripts and their
  // translations to account for the re-mapping process

  for (i=0;i<Gene_getTranscriptCount(gene);i++) {
    Transcript *trans = Gene_getTranscriptAt(gene,i);

    // need to grab the translation before starting to
    // re-jiggle the exons

// Transcript transforming is ODD
    Transcript_transform(trans, exonTransforms);
  }

  // unset the start, end, and strand - they need to be recalculated

#ifdef DONE
  $self->{_chr_name} = undef;
#endif

  IDHash_free(exonTransforms, NULL);

  return gene;
}

void Gene_free(Gene *gene) {
// NIY
  Object_decRefCount(gene);

  if (Object_getRefCount(gene) > 0) {
    return;
  } else if (Object_getRefCount(gene) < 0) {
    fprintf(stderr,"Error: Negative reference count for Gene\n"
                   "       Freeing it anyway\n");
  }
  printf("Gene_free not implemented\n");
}
