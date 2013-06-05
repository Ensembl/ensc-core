#define __GENE_MAIN__
#include "Gene.h"
#undef __GENE_MAIN__

#include "DBAdaptor.h"
#include "GeneAdaptor.h"
#include "AttributeAdaptor.h"
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

  gene->transcripts = Vector_new();

  gene->objectType = CLASS_GENE;
  Object_incRefCount(gene);

  gene->funcs = &geneFuncs;

  return gene;
}

Gene *Gene_shallowCopy(Gene *gene) {
  Gene *newGene = Gene_new();

  memcpy(newGene,gene,sizeof(Gene));

  return newGene;
}

/*
=head2 get_all_Attributes

  Arg [1]    : (optional) String $attrib_code
               The code of the attribute type to retrieve values for
  Example    : my ($author) = @{ $gene->get_all_Attributes('author') };
               my @gene_attributes = @{ $gene->get_all_Attributes };
  Description: Gets a list of Attributes of this gene.
               Optionally just get Attributes for given code.
  Returntype : Listref of Bio::EnsEMBL::Attribute
  Exceptions : warning if gene does not have attached adaptor and attempts lazy
               load.
  Caller     : general
  Status     : Stable

=cut
*/
// New

// NIY:
// Because this can filter the results the vector that gets returned must be freeable - so for now
// make a copy of the translation->attributes vector if returning unfiltered so behaviour is
// consistent. Long term probably want reference count incremented
Vector *Gene_getAllAttributes(Gene *gene, char *attribCode) {
  if (gene->attributes == NULL) {
    GeneAdaptor *ga = (GeneAdaptor *)Gene_getAdaptor(gene);
    if (ga == NULL) { // No adaptor
// Perl comments out the warning, I'll put it back for now, just in case
      fprintf(stderr,"Warning: Cannot get attributes without an adaptor.\n");
      return Vector_new();
    }

    AttributeAdaptor *ata = DBAdaptor_getAttributeAdaptor(ga->dba);
    gene->attributes = AttributeAdaptor_fetchAllByGene(ata, gene, NULL);
  }

  if (attribCode != NULL) {
    Vector *results = Vector_new();
    int i;
    for (i=0; i<Vector_getNumElement(gene->attributes); i++) {
      Attribute *attrib = Vector_getElementAt(gene->attributes, i);
      if (!strcasecmp(attrib->code, attribCode)) {
        Vector_addElement(results, attrib);
      }
    }
    return results;
  } else {
// See NIY note above for why I'm making a copy
    return Vector_copy(gene->attributes);
  }
}

// Note this used to be called getAllDBLinks but it looks more like the current getAllDBEntries so call it that
Vector *Gene_getAllDBEntries(Gene *g) {
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
//  Vector *exonVector = Vector_new();
//  void **values;

  for (i=0;i<Gene_getTranscriptCount(gene);i++) {
    Transcript *trans = Gene_getTranscriptAt(gene,i);
    int j;

    for (j=0;j<Transcript_getExonCount(trans);j++) {
      Exon *exon = Transcript_getExonAt(trans,j);
      if (!IDHash_contains(exonHash,(IDType)exon)) {
        IDHash_add(exonHash,(IDType)exon,exon);
      }
    }
  }

//  values = IDHash_getValues(exonHash);
//
//  for (i=0;i<IDHash_getNumValues(exonHash);i++) {
//    Vector_addElement(exonVector, values[i]);
//  }
//  free(values);
  Vector *exonVector = IDHash_getValuesVector(exonHash);
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
//  printf("Gene_free not implemented\n");
}
