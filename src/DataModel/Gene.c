#define __GENE_MAIN__
#include "Gene.h"
#undef __GENE_MAIN__
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

int Gene_setStart(Gene *gene, int start) {
  gene->start = start;
  Gene_setStartIsSet(gene,TRUE);
  return gene->start;
}

int Gene_getStart(Gene *gene) {
  int multiFlag = 0;

  if (Gene_getStartIsSet(gene) == FALSE) {
    char *lastContig = NULL;
    int i;
    Vector *exonVector = Gene_getAllExons(gene);

    for (i=0; i<Vector_getNumElement(exonVector); i++) {
      Exon *exon = Vector_getElementAt(exonVector,i);

      if (!Gene_getStartIsSet(gene) || 
          Exon_getStart(exon) < gene->start) {
        gene->start = Exon_getStart(exon);
        Gene_setStartIsSet(gene,TRUE);
      }
      if (multiFlag || 
          (lastContig && strcmp(lastContig, BaseContig_getName(Exon_getContig(exon))))) {
        multiFlag = 1;
      }
      lastContig =  BaseContig_getName(Exon_getContig(exon));
    }
    Vector_free(exonVector,NULL);
  }

  if (multiFlag) {
    fprintf(stderr, "WARNING: Gene_getStart - Gene spans multiple contigs."
                "The return value from getStart may not be what you want");
  }

  return gene->start;
}

int Gene_setEnd(Gene *gene, int end) {
  gene->end = end;
  Gene_setEndIsSet(gene,TRUE);
  return gene->end;
}

int Gene_getEnd(Gene *gene) {
  int multiFlag = 0;

  if (Gene_getEndIsSet(gene) == FALSE) {
    char *lastContig = NULL;
    int i;
    Vector *exonVector = Gene_getAllExons(gene);

    for (i=0; i<Vector_getNumElement(exonVector); i++) {
      Exon *exon = Vector_getElementAt(exonVector,i);

      if (!Gene_getEndIsSet(gene) || 
          Exon_getEnd(exon) > gene->end) {
        gene->end = Exon_getEnd(exon);
        Gene_setEndIsSet(gene,TRUE);
      }
      if (multiFlag || 
          (lastContig && strcmp(lastContig, BaseContig_getName(Exon_getContig(exon))))) {
        multiFlag = 1;
      }
      lastContig =  BaseContig_getName(Exon_getContig(exon));
    }
    Vector_free(exonVector,NULL);
  }

  if (multiFlag) {
    fprintf(stderr, "WARNING: Gene_getEnd - Gene spans multiple contigs."
                "The return value from getEnd may not be what you want");
  }

  return gene->end;
}

int Gene_setStrand(Gene *gene, int strand) {
  gene->strand = strand;
  Gene_setStrandIsSet(gene,TRUE);
  return gene->strand;
}

int Gene_getStrand(Gene *gene) {
  int multiFlag = 0;

  if (Gene_getStrandIsSet(gene) == FALSE) {
    char *lastContig = NULL;
    int i;
    Vector *exonVector = Gene_getAllExons(gene);

    if (!Vector_getNumElement(exonVector)) {
      fprintf(stderr, "WARNING: Gene_getStrand - Gene contains no exons."
                  "The return value from getStrand may not be what you want");
      return 0;
    }

    for (i=0; i<Vector_getNumElement(exonVector); i++) {
      Exon *exon = Vector_getElementAt(exonVector,i);

      if (lastContig && strcmp(lastContig, BaseContig_getName(Exon_getContig(exon)))) {
        multiFlag = 1;
        break;
      }
      lastContig =  BaseContig_getName(Exon_getContig(exon));
    }


    if (multiFlag) {
      fprintf(stderr, "WARNING: Gene_getStrand - Gene spans multiple contigs."
                  "The return value from getStrand may not be what you want");
      Vector_free(exonVector,NULL);
      return 0;
    }

    Gene_setStrand(gene, Exon_getStrand((Exon *)Vector_getElementAt(exonVector,0)));
    Vector_free(exonVector,NULL);
  }

  return gene->strand;
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
  Gene_setStartIsSet(gene,FALSE);
  Gene_setEndIsSet(gene,FALSE);
  Gene_setStrandIsSet(gene,FALSE);

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
  Gene_setStartIsSet(gene,FALSE);
  Gene_setEndIsSet(gene,FALSE);
  Gene_setStrandIsSet(gene,FALSE);

#ifdef DONE
  $self->{_chr_name} = undef;
#endif

  IDHash_free(exonTransforms, NULL);

  return gene;
}
