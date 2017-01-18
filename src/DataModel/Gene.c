/*
 * Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
 * Copyright [2016-2017] EMBL-European Bioinformatics Institute
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

  // Default gene biotype to protein_coding (what perl does)
  Gene_setBiotype(gene, "protein_coding");

  return gene;
}

Gene *Gene_shallowCopy(Gene *gene) {
  Gene *newGene = Gene_new();

  memcpy(newGene,gene,sizeof(Gene));

  return newGene;
}

/*
=head2 transform

  Arg [1]    : String - coordinate system name to transform to
  Arg [2]    : String - coordinate system version
  Example    : my $new_gene = $gene->transform('supercontig');
  Description: Moves this gene to the given coordinate system. If this gene has
               Transcripts attached, they move as well.
  Returntype : Bio::EnsEMBL::Gene
  Exceptions : throw on wrong parameters
  Caller     : general
  Status     : Stable

=cut
*/
// New NIY

Gene *Gene_transform(Gene *gene, char *csName, char *csVersion, Slice *toSlice) {
  fprintf(stderr, "Error: Gene_transform not implemented yet\n");
  exit(1);

/*
  my $self = shift;

  # catch for old style transform calls
  if( !@_  || ( ref $_[0] && ($_[0]->isa( "Bio::EnsEMBL::Slice" ) or $_[0]->isa( "Bio::EnsEMBL::LRGSlice" )) )) {
    deprecate('Calling transform without a coord system name is deprecated.');
    return $self->_deprecated_transform(@_);
  }

  my $new_gene = $self->SUPER::transform(@_);

  if ( !defined($new_gene) ) {
    # check if this gene projects at all to requested coord system,
    #  if not we are done.
    my @segments = @{ $self->project(@_) };
    if ( !@segments ) {
      return undef;
    }
  }

  #
  # If you are transforming the gene then make sure the transcripts and exons are loaded
  #

  foreach my $tran (@{$self->get_all_Transcripts}){
    $tran->get_all_Exons();
  }

  if( exists $self->{'_transcript_array'} ) {
    my @new_transcripts;
    my ( $strand, $slice );
    my $low_start = POSIX::INT_MAX;
    my $hi_end = POSIX::INT_MIN;
    for my $old_transcript ( @{$self->{'_transcript_array'}} ) {
      my $new_transcript = $old_transcript->transform( @_ );
      # this can fail if gene transform failed

      return undef unless $new_transcript;

      if( ! defined $new_gene ) {
        if( $new_transcript->start() < $low_start ) {
          $low_start = $new_transcript->start();
        }
        if( $new_transcript->end() > $hi_end ) {
          $hi_end = $new_transcript->end();
        }
        $slice = $new_transcript->slice();
        $strand = $new_transcript->strand();
      }
      push( @new_transcripts, $new_transcript );
    }

    if( ! defined $new_gene ) {
      %$new_gene = %$self;
      bless $new_gene, ref( $self );

      $new_gene->start( $low_start );
      $new_gene->end( $hi_end );
      $new_gene->strand( $strand );
      $new_gene->slice( $slice );
    }

    $new_gene->{'_transcript_array'} = \@new_transcripts;
  }

  if(exists $self->{attributes}) {
    $new_gene->{attributes} = [@{$self->{attributes}}];
  }

  return $new_gene;
*/
}



/*
=head2 transfer

  Arg [1]    : Bio::EnsEMBL::Slice $destination_slice
  Example    : my $new_gene = $gene->transfer($slice);
  Description: Moves this Gene to given target slice coordinates. If Transcripts
               are attached they are moved as well. Returns a new gene.
  Returntype : Bio::EnsEMBL::Gene
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/
// New
Gene *Gene_transfer(Gene *gene, Slice *slice) {
  // Call super transfer
  Gene *newGene = (Gene*)SeqFeature_transfer((SeqFeature*)gene, slice);

  if (newGene == NULL) {
    return NULL;
  }

  if (gene->transcripts &&  Vector_getNumElement(gene->transcripts)) {
    Vector *newTranscripts = Vector_new();
    int i;
    for (i=0; i<Vector_getNumElement(gene->transcripts); i++) {
      Transcript *oldTranscript = Vector_getElementAt(gene->transcripts, i);
      Transcript *newTranscript = Transcript_transfer(oldTranscript, slice);
      Vector_addElement(newTranscripts, newTranscript);
    }
    newGene->transcripts = newTranscripts;
  }

/* Seems very hacky - don't do for now
  if (exists $self->{attributes}) {
    $new_gene->{attributes} = [@{$self->{attributes}}];
  }
*/

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
//      fprintf(stderr,"Warning: Cannot get attributes without an adaptor.\n");
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
  fprintf(stderr, "DEPRECATED: This column has been removed in 74.\n  It does not set anything and it just returns NULL\n");
  return NULL;
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

int Gene_getExonCount(Gene *gene) {
  IDHash *exonHash = IDHash_new(IDHASH_SMALL);
  int i;
  int nExon = 0;

  for (i=0;i<Gene_getTranscriptCount(gene);i++) {
    Transcript *trans = Gene_getTranscriptAt(gene,i);
    int j;

    for (j=0;j<Transcript_getExonCount(trans);j++) {
      Exon *exon = Transcript_getExonAt(trans,j);
      if (!IDHash_contains(exonHash,(IDType)exon)) {
        IDHash_add(exonHash,(IDType)exon,exon);
        nExon++;
      }
    }
  }

  IDHash_free(exonHash,NULL);
  
  return nExon;
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
//  fprintf(stderr,"Gene_free called\n");

  if (Object_getRefCount(gene) > 0) {
    //fprintf(stderr,"Not freeing gene with ref count = %d\n", Object_getRefCount(gene));
    return;
  } else if (Object_getRefCount(gene) < 0) {
    fprintf(stderr,"Error: Negative reference count for Gene\n"
                   "       Freeing it anyway\n");
  }

  //fprintf(stderr,"  Freeing gene %p %s\n", gene, Gene_getStableId(gene));
  int i;
  for (i=0; i<Gene_getTranscriptCount(gene); i++) {
    Transcript *trans = Gene_getTranscriptAt(gene, i);
    Transcript_free(trans);
  }

  Vector_free(gene->transcripts);

  StableIdInfo_freePtrs(&gene->si);

  free(gene);
}

/*
=head2 add_Transcript

  Arg [1]    : Bio::EnsEMBL::Transcript $trans
               The transcript to add to the gene
  Example    : my $transcript = Bio::EnsEMBL::Transcript->new(...);
               $gene->add_Transcript($transcript);
  Description: Adds another Transcript to the set of alternatively
               spliced Transcripts of this gene. If it shares exons
               with another Transcript, these should be object-identical.
  Returntype : none
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/
// New
void Gene_addTranscript(Gene *gene, Transcript *trans) {

  if (gene->transcripts == NULL) {
    gene->transcripts = Vector_new();
  }

  Vector_addElement(gene->transcripts, trans);

  Gene_recalculateCoordinates(gene);
}


/*
=head2 recalculate_coordinates

  Example    : $gene->recalculate_coordinates;
  Description: Called when transcript added to the gene, tries to adapt the
               coords for the gene.
  Returntype : none
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut
*/
// New
void Gene_recalculateCoordinates(Gene *gene) {
  Vector *transcripts = Gene_getAllTranscripts(gene);
  
  if (transcripts == NULL || Vector_getNumElement(transcripts) == 0) {
    return;
  }

  Slice *slice;
  long start;
  long end;
  int strand;

  Transcript *firstTrans = Vector_getElementAt(transcripts, 0);

  slice  = Transcript_getSlice(firstTrans);
  strand = Transcript_getStrand(firstTrans);
  start  = Transcript_getStart(firstTrans);
  end    = Transcript_getEnd(firstTrans);

  int transSplicing = 0;

  int i;
  for (i=0; i<Vector_getNumElement(transcripts); i++) {
    Transcript *t = Vector_getElementAt(transcripts, i);
    if (Transcript_getStart(t) < start) {
      start = Transcript_getStart(t);
    }

    if (Transcript_getEnd(t) > end) {
      end = Transcript_getEnd(t);
    }

    if (EcoString_strcmp(Slice_getName(slice), Slice_getName(Transcript_getSlice(t)))) {
      fprintf(stderr, "Transcripts with different slices not allowed on one Gene\n");
      exit(1);
    }

    if ( Transcript_getStrand(t) != strand ) {
      transSplicing = 1;
    }
  }

  if (transSplicing) {
    fprintf(stderr,"Warning: Gene contained trans splicing event\n");
  }

  Gene_setStart(gene, start);
  Gene_setEnd(gene, end);
  Gene_setStrand(gene, strand);
  Gene_setSlice(gene, slice);
}


/*
=head2 get_all_Transcripts

  Example    : my @transcripts = @{ $gene->get_all_Transcripts };
  Description: Returns the Transcripts in this gene.
  Returntype : Listref of Bio::EnsEMBL::Transcript objects
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/
// New
// NIY Transcript lazy loading
Vector *Gene_getAllTranscripts(Gene *gene) {

/*
  if( ! exists $self->{'_transcript_array'} ) {
    if( defined $self->adaptor() ) {
      my $ta = $self->adaptor()->db()->get_TranscriptAdaptor();
      my $transcripts = $ta->fetch_all_by_Gene( $self );
      $self->{'_transcript_array'} = $transcripts;
    }
  }
*/
  return gene->transcripts;
}






