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
  Gene_setIsCurrent(gene,1);

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
//    GeneAdaptor_getStableEntryInfo(ga,gene);
    fprintf(stderr, "New Gene code shouldn't need to lazy load stable ids\n");
    exit(1);
  }
  return StableIdInfo_getStableId(&(gene->si));
}

char *Gene_setDescription(Gene *g, char *description) {
  if ((g->description = (char *)malloc(strlen(description)+1)) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for description\n");
    return NULL;
  }

  strcpy(g->description,description);

  return g->description;
}

char *Gene_setCanonicalAnnotation(Gene *g, char *canonicalAnnotation) {
  if ((g->canonicalAnnotation = (char *)malloc(strlen(canonicalAnnotation)+1)) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for canonicalAnnotation\n");
    return NULL;
  }

  strcpy(g->canonicalAnnotation,canonicalAnnotation);

  return g->canonicalAnnotation;
}


ECOSTRING Gene_setBiotype(Gene *g, char *biotype) {
  EcoString_copyStr(ecoSTable, &(g->biotype),biotype,0);

  if (g->biotype == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for biotype\n");
    return NULL;
  }

  return g->biotype;
}

ECOSTRING Gene_setStatus(Gene *g, char *status) {
  EcoString_copyStr(ecoSTable, &(g->status),status,0);

  if (g->status == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for status\n");
    return NULL;
  }

  return g->status;
}

ECOSTRING Gene_setSource(Gene *g, char *source) {
  EcoString_copyStr(ecoSTable, &(g->source),source,0);

  if (g->source == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for source\n");
    return NULL;
  }

  return g->source;
}

ECOSTRING Gene_setExternalDb(Gene *g, char *externalDb) {
  EcoString_copyStr(ecoSTable, &(g->externalDb),externalDb,0);

  if (g->externalDb == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for externalDb\n");
    return NULL;
  }

  return g->externalDb;
}

ECOSTRING Gene_setExternalStatus(Gene *g, char *externalStatus) {
  EcoString_copyStr(ecoSTable, &(g->externalStatus),externalStatus,0);

  if (g->externalStatus == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for externalStatus\n");
    return NULL;
  }

  return g->externalStatus;
}

char *Gene_setExternalName(Gene *g, char *externalName) {
  if (externalName == NULL) {
    g->externalName = NULL;
    return NULL;
  }
  if ((g->externalName = (char *)malloc(strlen(externalName)+1)) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for externalName\n");
    return NULL;
  }

  strcpy(g->externalName,externalName);

  return g->externalName;
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
