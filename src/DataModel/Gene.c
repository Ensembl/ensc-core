#include "Gene.h"
#include "GeneAdaptor.h"
#include "IDHash.h"

Gene *Gene_new() {
  Gene *gene;

  if ((gene = (Gene *)calloc(1,sizeof(Gene))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for gene\n");
    return NULL;
  }

  return gene;
}

char *Gene_getStableId(Gene *gene) {
  GeneAdaptor *ga = (GeneAdaptor *)Gene_getAdaptor(gene);

  if (StableIdInfo_getStableId(&(gene->si)) == NULL && ga) {
    GeneAdaptor_getStableEntryInfo(ga,gene);
  }
  return StableIdInfo_getStableId(&(gene->si));
}

char *Gene_setType(Gene *g, char *type) {
  if ((g->type = (char *)malloc(strlen(type))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for gene type\n");
    return NULL;
  }

  strcpy(g->type,type);

  return g->type;
}

void Gene_transformToSlice(Gene *gene, Slice *slice) {
  IDHash *exonTransforms = IDHash_new(IDHASH_SMALL);

#ifdef DONE
  Exon = Gene_getAllExons(gene);

  // transform Exons
  for my $exon ( @{$self->get_all_Exons()} ) {
    my $newExon = $exon->transform( $slice );
    $exon_transforms{ $exon } = $newExon;
  }

  // now need to re-jiggle the transcripts and their
  // translations to account for the re-mapping process

  foreach my $transcript ( @{$self->get_all_Transcripts()} ) {

    // need to grab the translation before starting to
    // re-jiggle the exons

    $transcript->transform( \%exon_transforms );
  }

  // unset the start, end, and strand - they need to be recalculated
  $self->{start} = undef;
  $self->{end} = undef;
  $self->{strand} = undef;
  $self->{_chr_name} = undef;
#endif

  IDHash_free(exonTransforms, NULL);

  return;
}
