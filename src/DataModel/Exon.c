/*
 * Copyright [1999-2017] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

#define __EXON_MAIN__
#include "Exon.h"
#undef __EXON_MAIN__
#include "ExonAdaptor.h"
#include "AssemblyMapperAdaptor.h"
#include "SupportingFeatureAdaptor.h"
#include "AssemblyMapper.h"
#include "DBAdaptor.h"
#include "SliceAdaptor.h"
#include "ChromosomeAdaptor.h"
#include "RawContigAdaptor.h"
#include "StrUtil.h"
#include "CoordSystem.h"
#include "SeqUtil.h"
#include "CachingSequenceAdaptor.h"
#include "DNAAlignFeature.h"
#include "BaseAlignFeature.h"
#include "Object.h"
#include "translate.h"
#include "StableIdInfo.h"

Exon *Exon_new() {
  Exon *exon;

  if ((exon = (Exon *)calloc(1,sizeof(Exon))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for exon\n");
    return NULL;
  }

  exon->objectType = CLASS_EXON;
//  Object_incRefCount(exon);

  exon->funcs = &exonFuncs;

/* Set to empty values */
  Exon_setModified(exon,0);
  Exon_setCreated(exon,0);
  Exon_setVersion(exon,-1);
  Exon_setIsCurrent(exon,1);

  Exon_setStart(exon, POS_UNDEF);
  Exon_setEnd(exon, POS_UNDEF);
  Exon_setStrand(exon, STRAND_UNDEF);
  Exon_setPhase(exon, PHASE_UNDEF);
  Exon_setEndPhase(exon, PHASE_UNDEF);


  return exon;
}

/*
=head2 hashkey

  Arg [1]    : none
  Example    : if(exists $hash{$exon->hashkey}) { do_something(); }
  Description: Returns a unique hashkey that can be used to uniquely identify
               this exon.  Exons are considered to be identical if they share
               the same seq_region, start, end, strand, phase, end_phase.
               Note that this will consider two exons on different slices
               to be different, even if they actually are not.
  Returntype : string formatted as slice_name-start-end-strand-phase-end_phase
  Exceptions : thrown if not all the necessary attributes needed to generate
               a unique hash value are set
               set
  Caller     : general
  Status     : Stable

=cut
*/
// New
// Adapted to take empty hashKey string as arg so don't have to allocate
void Exon_getHashKey(Exon *exon, char *hashKey) {
  Slice *slice    = Exon_getSlice(exon);
  char *sliceName = slice != NULL ? Slice_getName(slice) : NULL;
  long start      = Exon_getStart(exon);
  long end        = Exon_getEnd(exon);
  int strand      = Exon_getStrand(exon);
  int phase       = Exon_getPhase(exon);
  int endPhase    = Exon_getEndPhase(exon);


  if (sliceName == NULL) {
    fprintf(stderr, "Slice must be set to generate correct exon hashkey.\n");
    exit(1);
  }

  if (start == POS_UNDEF) {
    fprintf(stderr, "start attribute must be defined to generate correct hashkey.");
  }

  if (end == POS_UNDEF) {
    fprintf(stderr, "end attribute must be defined to generate correct hashkey.");
  }

  if (strand == STRAND_UNDEF) {
    fprintf(stderr, "strand attribute must be defined to generate correct hashkey.");
  }

  if (phase == PHASE_UNDEF) {
    fprintf(stderr, "phase attribute must be defined to generate correct hashkey.");
  }

  if (endPhase == PHASE_UNDEF) {
    fprintf(stderr, "endPhase attribute must be defined to generate correct hashkey.");
  }

  sprintf(hashKey, "%s-%ld-%ld-%d-%d-%d", sliceName, start, end, strand, phase, endPhase);

  return;
}

/*
=head2 flush_supporting_features

  Example     : $exon->flush_supporting_features;
  Description : Removes all supporting evidence from the exon.
  Return type : (Empty) listref
  Exceptions  : none
  Caller      : general
  Status      : Stable

=cut
*/
// New
void Exon_flushSupportingFeatures(Exon *exon) {
  if (exon->supportingFeatures != NULL) {
// NIY: Do we need to free features??
    //exon->supportingFeatures->freeElement = NULL;
    Vector_setFreeFunc(exon->supportingFeatures, Object_freeImpl);
    Vector_free(exon->supportingFeatures);
  }

  exon->supportingFeatures = Vector_new();
}


Exon *Exon_shallowCopyImpl(Exon *exon) {
  Exon *newExon = Exon_new();

  memcpy(newExon,exon,sizeof(Exon));

// Make a copy of this string
  
  if (exon->si.stableId) StrUtil_copyString(&(newExon->si.stableId), exon->si.stableId, 0);

  return newExon;
}

char *Exon_getStableId(Exon *exon) {
  ExonAdaptor *ea = (ExonAdaptor *)Exon_getAdaptor(exon);

  if (StableIdInfo_getStableId(&(exon->si)) == NULL && ea) {
//    ExonAdaptor_getStableEntryInfo(ea,exon);
  }
  return StableIdInfo_getStableId(&(exon->si));
}

time_t Exon_getCreated(Exon *exon) {
  ExonAdaptor *ea = (ExonAdaptor *)Exon_getAdaptor(exon);

  if (StableIdInfo_getCreated(&(exon->si)) == 0 && ea) {
//    ExonAdaptor_getStableEntryInfo(ea,exon);
  }
  return StableIdInfo_getCreated(&(exon->si));
}

time_t Exon_getModified(Exon *exon) {
  ExonAdaptor *ea = (ExonAdaptor *)Exon_getAdaptor(exon);

  if (StableIdInfo_getModified(&(exon->si)) == 0 && ea) {
//    ExonAdaptor_getStableEntryInfo(ea,exon);
  }
  return StableIdInfo_getModified(&(exon->si));
}

int Exon_getVersion(Exon *exon) {
  ExonAdaptor *ea = (ExonAdaptor *)Exon_getAdaptor(exon);

  if (StableIdInfo_getVersion(&(exon->si)) == -1 && ea) {
//    ExonAdaptor_getStableEntryInfo(ea,exon);
  }
  return StableIdInfo_getVersion(&(exon->si));
}

int Exon_forwardStrandCompFunc(const void *a, const void *b) {
  Exon **e1 = (Exon **)a;
  Exon **e2 = (Exon **)b;

  if (Exon_getStart(*e1) > Exon_getStart(*e2)) {
    return 1;
  } else if (Exon_getStart(*e1) < Exon_getStart(*e2)) {
    return -1;
  } else {
    return 0;
  }
}

int Exon_reverseStrandCompFunc(const void *a, const void *b) {
  Exon *e1 = *((Exon **)a);
  Exon *e2 = *((Exon **)b);

  return Exon_getStart(e2) - Exon_getStart(e1);
}

/*
=head2 transfer

  Arg [1]    : Bio::EnsEMBL::Slice $destination_slice
  Example    : none
  Description: Moves this Exon to given target slice coordinates. If Features
               are attached they are moved as well. Returns a new exon.
  Returntype : Bio::EnsEMBL::Gene
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/
// New
static int notWarned = 1;
Exon *Exon_transfer(Exon *exon, Slice *slice) {

  // Call super transfer
  //my $new_exon = $self->SUPER::transfer( @_ );
  Exon *newExon = (Exon*)SeqFeature_transfer((SeqFeature*)exon, slice);

  if (newExon == NULL) {
    return NULL;
  }

  if (notWarned) {
    fprintf(stderr,"seq cache clearing not implemented yet for exon\n");
    notWarned = 0;
  }
  if (exon->supportingFeatures != NULL && Vector_getNumElement(exon->supportingFeatures) != 0) {
    Vector *newFeatures = Vector_new();

    int i;
    for (i=0; i<Vector_getNumElement(exon->supportingFeatures); i++) {
      SeqFeature *oldFeature = Vector_getElementAt(exon->supportingFeatures, i);
      SeqFeature *newFeature = SeqFeature_transfer(oldFeature, slice);
      Vector_addElement(newFeatures, newFeature);
    }
    newExon->supportingFeatures = newFeatures;
  }

/*
  #dont want to share the same sequence cache
  delete $new_exon->{'_seq_cache'};
*/
  newExon->seqCacheString = NULL;

  return newExon;
}



/*
=head2 coding_region_start

    Arg [1]     : Bio::EnsEMBL::Transcript $transcript
    Example     : $coding_region_start =
                    $exon->coding_region_start($transcript);
    Description : Returns the start position of the coding region
                  of the exon in slice-relative coordinates on the
                  forward strand.  Returns undef if the whole exon is
                  non-coding.
                  Since an exon may be part of one or more transcripts,
                  the relevant transcript must be given as argument to
                  this method.
    Return type : Integer or undef
    Exceptions  : Throws if the given argument is not a transcript.
    Caller      : General
    Status      : Stable

=cut

# The implementation of this method is analogous to the implementation
# of cdna_coding_start().
*/
// New
long Exon_getCodingRegionStart(Exon *exon, Transcript *transcript) {

/* NIY
  if (defined($id) && /exists $self->{coding_region_start}->{$id}) {
    return $self->{coding_region_start}->{$id};
  }
*/
  //fprintf(stderr, "Caching of Exon_codingRegionStart not implemented yet\n");

  long codingRegionStart;
  long transcriptCodingStart = Transcript_getCodingRegionStart(transcript);
  if (transcriptCodingStart != POS_UNDEF) {
    long start = Exon_getStart(exon);

    if (transcriptCodingStart < start) {
      // Coding region starts upstream of this exon...

      if (Transcript_getCodingRegionEnd(transcript) < start) {
        // ... and also ends upstream of this exon.
        codingRegionStart = POS_UNDEF;
      } else {
        // ... and does not end upstream of this exon.
        codingRegionStart = start;
      }
    } else {
      // Coding region starts either within or downstream of this
      // exon.

      if (transcriptCodingStart <= Exon_getEnd(exon)) {
        // Coding region starts within this exon.
        codingRegionStart = transcriptCodingStart;
      } else {
        // Coding region starts downstream of this exon.
        codingRegionStart = POS_UNDEF;
      }
    }
  } else {
    codingRegionStart = POS_UNDEF;
  }

/* NIY
  if (defined $id) {
    $self->{codingRegionStart}->{$id} = $codingRegionStart;
    $self->{coding_region_end}->{$id} = undef if ! defined $codingRegionStart;
  }
*/

  return codingRegionStart;
}

/*
=head2 coding_region_end

    Arg [1]     : Bio::EnsEMBL::Transcript $transcript
    Example     : $coding_region_end =
                    $exon->coding_region_end($transcript);
    Description : Returns the end position of the coding region of
                  the exon in slice-relative coordinates on the
                  forward strand.  Returns undef if the whole exon is
                  non-coding.
                  Since an exon may be part of one or more transcripts,
                  the relevant transcript must be given as argument to
                  this method.
    Return type : Integer or undef
    Exceptions  : Throws if the given argument is not a transcript.
    Caller      : General
    Status      : Stable

=cut

# The implementation of this method is analogous to the implementation
# of cdna_coding_end().
*/
// New
long Exon_getCodingRegionEnd(Exon *exon, Transcript *transcript) {

/* NIY
  if(defined $id && exists $self->{coding_region_end}->{$id}) {
    return $self->{coding_region_end}->{$id};
  }
*/
  //fprintf(stderr, "Caching of Exon_codingRegionEnd not implemented yet\n");

  long codingRegionEnd;
  long transcriptCodingEnd = Transcript_getCodingRegionEnd(transcript);
  if (transcriptCodingEnd != POS_UNDEF) {

    long end = Exon_getEnd(exon);
    if (transcriptCodingEnd > end) {
      // Coding region ends downstream of this exon...

      if (Transcript_getCodingRegionStart(transcript) > end) {
        // ... and also starts downstream of this exon.
        codingRegionEnd = POS_UNDEF;
      } else {
        // ... and does not start downstream of this exon.
        codingRegionEnd = end;
      }
    } else {
      // Coding region ends either within or upstream of this
      // exon.
      if (transcriptCodingEnd >= Exon_getStart(exon)) {
        codingRegionEnd = transcriptCodingEnd;
      } else {
        codingRegionEnd = POS_UNDEF;
      }
    }
  } else {
    // This is a non-coding transcript.
    codingRegionEnd = POS_UNDEF;
  }

/* NIY
  if(defined $id) {
    $self->{codingRegionEnd}->{$id} = $codingRegionEnd;
    $self->{coding_region_start}->{$id} = undef if ! defined $codingRegionEnd;
  }
*/

  return codingRegionEnd;
}



Exon *Exon_transformRawContigToSliceImpl(Exon *exon, Slice *slice) {
  Exon *newExon;
  BaseAdaptor *adaptor;
  AssemblyMapperAdaptor *ama;
  AssemblyMapper *assMapper;
  MapperRangeSet *mapped;
  MapperCoordinate *mc;

  if (!Exon_getContig(exon)) {
    fprintf(stderr,"ERROR: Exon's contig must be defined to transform to Slice coords");
    exit(1);
  }
  // fprintf(stderr,"transforming exon %d from raw contig to slice coords\n",Exon_getDbID(exon));
  // fprintf(stderr,"exon %s\n",Exon_getStableId(exon));

  adaptor = Slice_getAdaptor(slice);
  if (!adaptor) {
    adaptor = BaseContig_getAdaptor((BaseContig *)Exon_getContig(exon));
  }
  
  if (!adaptor) {
    fprintf(stderr, "ERROR: Cannot transform to exon slice unless either the " 
		    "exon->contig->adaptor or slice->adaptor is defined");
    exit(1);
  }

  ama = DBAdaptor_getAssemblyMapperAdaptor(adaptor->dba);
  assMapper = AssemblyMapperAdaptor_fetchByType(ama, CoordSystem_getName(Slice_getCoordSystem(slice)));
  
  mapped = AssemblyMapper_mapCoordinatesToAssembly(assMapper,
     BaseContig_getName(Exon_getContig(exon)),//BaseContig_getDbID(Exon_getContig(exon)),
     Exon_getStart(exon),
     Exon_getEnd(exon),
     Exon_getStrand(exon)
    );

  // exons should always transform so in theory no error check necessary
  // actually we could have exons inside and outside the Slice 
  // because of db design and the query that produces them
  if (!mapped || mapped->nRange == 0) {
    fprintf(stderr, "ERROR: Exon couldnt map dbID = " IDFMTSTR "\n",Exon_getDbID(exon));
    exit(1);
  }

  // should get a gap object returned if an exon lies outside of 
  // the current slice.  Simply return the exon as is - i.e. untransformed.
  // this untransformed exon will be distinguishable as it will still have
  // contig attached to it and not a slice.
  if( MapperRangeSet_getRangeAt(mapped,0)->rangeType == MAPPERRANGE_GAP) {
    fprintf(stderr,"Exon in gap\n");
    return exon;
  }

  mc = (MapperCoordinate *)MapperRangeSet_getRangeAt(mapped,0);

  // fprintf(stderr," Exon %s mapped range = %d %d\n",Exon_getStableId(exon),mc->start,mc->end);

  // the slice is an empty slice, create an enitre chromosome slice and
  // replace the empty slice with it
//  if (Slice_getEmptyFlag(slice)) {
//    SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(adaptor->dba);
//    ChromosomeAdaptor *ca = DBAdaptor_getChromosomeAdaptor(adaptor->dba);
//    char *chrName = Chromosome_getName(ChromosomeAdaptor_fetchByDbID(ca,mc->id));

//    slice = SliceAdaptor_fetchByChrName(sa, chrName);
// NIY free old slice (or have special empty one)???
//  } 

  newExon = Exon_new();

  newExon = Exon_copy(newExon, exon, SHALLOW_DEPTH);

  if (Slice_getStrand(slice) == 1) {
    Exon_setStart(newExon,mc->start - Slice_getChrStart(slice) + 1);
    Exon_setEnd(newExon,mc->end - Slice_getChrStart(slice) + 1);
  } else {
    Exon_setStart(newExon,Slice_getChrEnd(slice) - mc->end + 1);
    Exon_setEnd(newExon,Slice_getChrEnd(slice) - mc->start + 1);
  }

  Exon_setStrand(newExon, mc->strand * Slice_getStrand(slice));
  Exon_setContig(newExon, slice);

  //copy the attached supporting features and transform them
#ifdef DONE
// NIY
  my @feats;
  if( exists $self->{_supporting_evidence} ) {
    foreach my $sf (@{$self->get_all_supporting_features()}) {
      push @feats, $sf->transform($slice);
    }
    $newexon->add_supporting_features(@feats);
  }
#endif

  // NIY free old exon (reference counting needed ???)

  // fprintf(stderr,"transformed exon %s %d %d\n",
  //        Exon_getStableId(newExon),
  //        Exon_getStart(newExon),
  //        Exon_getEnd(newExon)
  //       );
  return newExon;
}

Exon *Exon_copy(Exon *copy, Exon *orig, CopyDepth depth) {
  if (depth != SHALLOW_DEPTH) {
    fprintf(stderr, "ERROR: Only SHALLOW copy implemented in Exon\n");
    exit(1);
  }
  Exon_setStart(copy,Exon_getStart(orig));  
  Exon_setEnd(copy,Exon_getEnd(orig));  
  Exon_setDbID(copy,Exon_getDbID(orig));  
  Exon_setAdaptor(copy,Exon_getAdaptor(orig));  
  Exon_setStrand(copy,Exon_getStrand(orig));  
  Exon_setPhase(copy,Exon_getPhase(orig));  
  Exon_setEndPhase(copy,Exon_getEndPhase(orig));  
  Exon_setAnalysis(copy,Exon_getAnalysis(orig));  
  Exon_setContig(copy,Exon_getContig(orig));  

  return copy;
}

void Exon_findSupportingEvidence(Exon *exon, Vector *features, int isSorted) {
  int i;
  Vector *support = Vector_new();

  for (i=0; i<Vector_getNumElement(features); i++) {
    SeqFeature *f = Vector_getElementAt(features,i);

    // return if we have a sorted feature array
    if (isSorted == 1 && SeqFeature_getStart(f) > Exon_getEnd(exon)) {
      return;
    }
/* NIY subfeatures
    if ($f->sub_SeqFeature) {
      my @subf = $f->sub_SeqFeature;

      $self->find_supporting_evidence(\@subf);
    } else {
*/
/* NIY Speedup */
      if (!strcmp(BaseContig_getName(SeqFeature_getContig(f)), 
                  BaseContig_getName(Exon_getContig(exon)))) {
        if (SeqFeature_getEnd(f) >= Exon_getStart(exon) && 
            SeqFeature_getStart(f) <= Exon_getEnd(exon) && 
            SeqFeature_getStrand(f) == Exon_getStrand(exon)) {
          Vector_addElement(support,f);
        }
      }
/* NIY subfeatures
    }
*/
  }
  Exon_addSupportingFeatures(exon, support);

  Vector_free(support);
}

Vector *Exon_getAllSupportingFeaturesImpl(Exon *exon) {

  if (!exon->supportingFeatures)  {
    if (Exon_getAdaptor(exon)) {
      DBAdaptor *dba = Exon_getAdaptor(exon)->dba;
      SupportingFeatureAdaptor *sfa = DBAdaptor_getSupportingFeatureAdaptor(dba);
      exon->supportingFeatures = SupportingFeatureAdaptor_fetchAllByExon(sfa,exon);
    }
  }

  if (!exon->supportingFeatures) {
    exon->supportingFeatures = emptyVector;
  }
  return exon->supportingFeatures;
}

void Exon_addSupportingFeaturesImpl(Exon *exon, Vector *v) {
  if (!exon->supportingFeatures || exon->supportingFeatures == emptyVector) {
    exon->supportingFeatures = Vector_new();
  }
  int i;
  for (i=0; i<Vector_getNumElement(v); i++) {
    Object *obj = Vector_getElementAt(v, i);
    Object_incRefCount(obj);
  }
  Vector_append(exon->supportingFeatures,v);
}

void Exon_addSupportingFeature(Exon *exon, SeqFeature *sf) {
  if (!exon->supportingFeatures || exon->supportingFeatures == emptyVector) {
    exon->supportingFeatures = Vector_new();
  }
  //fprintf(stderr,"Adding support "IDFMTSTR" to exon "IDFMTSTR" %p\n",SeqFeature_getDbID(sf),Exon_getDbID(exon),exon);
  Object_incRefCount(sf);
  Vector_addElement(exon->supportingFeatures,sf);
}

Exon *Exon_adjustStartEndImpl(Exon *exon, int startAdjust, int endAdjust) {

  Exon *newExon = Exon_new();

// Copy - NIY won't copy support
  Exon_copy(newExon, exon, SHALLOW_DEPTH);

  // invalidate the sequence cache
  // NIY delete $new_exon->{'_seq_cache'};
  newExon->seqCacheString = NULL;

/*
  fprintf(stderr,"adjustStartEndImpl with exon %ld-%ld:%d adjustStart %d adjustEnd %d\n", 
          Exon_getStart(exon), Exon_getEnd(exon), Exon_getStrand(exon), startAdjust, endAdjust);
*/

  if (Exon_getStrand(exon) == 1 ) {
    Exon_setStart(newExon, Exon_getStart(exon) + startAdjust );
    Exon_setEnd(newExon, Exon_getEnd(exon) + endAdjust );
  } else {
    Exon_setStart(newExon, Exon_getStart(exon) - endAdjust );
    Exon_setEnd(newExon, Exon_getEnd(exon) - startAdjust );
  }

// NIY Delete old exon

  Object_incRefCount(newExon);

  return newExon;
}

char *Exon_getPeptideImpl(Exon *exon, Transcript *trans) {
  char *peptide;
  MapperRangeSet *mapped;
  Vector *coords;
  int i;

  if (!trans) {
    fprintf(stderr, "Error: transcript arg non null in getPeptide\n");
    exit(1);
  }

  // convert exons coordinates to peptide coordinates
  mapped = Transcript_genomic2Pep(trans, Exon_getStart(exon), Exon_getEnd(exon),
                                  Exon_getStrand(exon),Exon_getContig(exon));

  coords = Vector_new();

  // filter out gaps
  for (i=0;i<mapped->nRange;i++) {
    MapperRange *mr = MapperRangeSet_getRangeAt(mapped,i);
    if (mr->rangeType == MAPPERRANGE_COORD) {
      Vector_addElement(coords,mr);
    }
  }

  // if this is UTR then the peptide will be empty string
  if (Vector_getNumElement(coords) > 1) {
    fprintf(stderr, "Error. Exon maps to multiple locations in peptide."
                    " Is this exon [%p] a member of this transcript [%p]?",
                    exon,trans);
    exit(1);

  } else if (Vector_getNumElement(coords) == 1) {
    MapperCoordinate *c = Vector_getElementAt(coords,0);
    int start,end;
    char *wholePeptide = Transcript_translate(trans);
    int lenPeptide = strlen(wholePeptide);

    // bioperl doesn't give back residues for incomplete codons
    // make sure we don't subseq too far...

    end = (c->end > lenPeptide) ? lenPeptide : c->end;
    start = (c->start < end) ? c->start : end;
// NIY free
// NIY check for off by one
    peptide = StrUtil_substr(wholePeptide,start,(end-start+1));
  }

  Vector_free(coords);
  MapperRangeSet_free(mapped);

  return peptide;
}

#ifdef DONE
void Exon_setSeq(Exon *exon, Sequence *seq) {
  fprintf(stderr, "Warning: Exon seq setting not supported\n");
  return;
}

sub Exon_getSeq(Exon *exon) {
  my $self = shift;
  my $arg = shift;

  if( defined $arg ) {
    $self->warn( "seq setting on Exon not supported currently" );
    $self->{'_seq_cache'} = $arg->seq();
  }

  if( defined $self->{'_seq_cache'} ) {
    return Bio::Seq->new(-seq=> $self->{'_seq_cache'});
  }

  my $seq;

  if ( ! defined $self->contig ) {
    $self->warn(" this exon doesn't have a contig you won't get a seq \n");
    return undef;
  }
  else {

    $seq = $self->contig()->subseq($self->start, $self->end);

    if($self->strand == -1){
      $seq =~ tr/ATGCatgc/TACGtacg/;
      $seq = reverse($seq);
    }

   }
  $self->{'_seq_cache'} = $seq;

  return Bio::Seq->new(-seq     => $self->{'_seq_cache'},
                       -id      => $self->stable_id,
                       -moltype => 'dna');
}
#endif

char  *Exon_getSeqStringImpl(Exon *exon) {
  char *seq;

  if (Exon_getSeqCacheString(exon)) {
    //fprintf(stderr,"Returning cached exon seq str\n");
    return Exon_getSeqCacheString(exon);
  }

  if (!Exon_getContig(exon)) {
    fprintf(stderr, "Warning: this exon %s doesn't have a contig you won't get a seq\n", Exon_getStableId(exon));
    return NULL;
  } else {

/*
    seq = BaseContig_getSubSeq(Exon_getContig(exon), 
                               Exon_getStart(exon), 
                               Exon_getEnd(exon),
                               1);
*/
    DBAdaptor *dba = Slice_getAdaptor(Exon_getSlice(exon))->dba;
    CachingSequenceAdaptor *csa = DBAdaptor_getCachingSequenceAdaptor(dba);
    seq = CachingSequenceAdaptor_fetchBySliceStartEndStrand(csa, 
                                                            Exon_getSlice(exon),
                                                            Exon_getStart(exon),
                                                            Exon_getEnd(exon),
                                                            1);

    //fprintf(stderr,"fetched exon seq str %s\n", seq);
    if (Exon_getStrand(exon) == -1){
//      SeqUtil_reverseComplement(seq,strlen(seq));
      
      if (Exon_getLength(exon) > 1000000) {
        char *tmpSeq;
        if ((tmpSeq = (char *)calloc(Exon_getLength(exon) + 3, sizeof(char))) == NULL) {
          fprintf(stderr, "Failed allocating temporary buffer for exon rev comp seq string\n");
          exit(1);
        }
        rev_comp(seq, tmpSeq, Exon_getLength(exon));
        strcpy(seq,tmpSeq);
        free(tmpSeq);
      } else {
        char *tmpSeq = NULL;

        if ((tmpSeq = (char *)calloc(1000002,sizeof(char))) == NULL) {
          fprintf(stderr,"Failed allocating tmpSeq\n");
          return NULL;
        }

        rev_comp(seq, tmpSeq, Exon_getLength(exon));
        strcpy(seq,tmpSeq);

        free(tmpSeq);
      }
    }
  }
  Exon_setSeqCacheString(exon, seq);
  return Exon_getSeqCacheString(exon);
}

void Exon_setSeqCacheString(Exon *exon, char *seq) {
  //fprintf(stderr,"Setting seq cache string for exon %p %s to %p\n", exon, Exon_getStableId(exon), seq);
  exon->seqCacheString = seq;
}

void Exon_clearSeqCacheString(Exon *exon) {
  if (exon->seqCacheString) free(exon->seqCacheString);

  exon->seqCacheString = NULL;
}

Exon *Exon_transformSliceToRawContigImpl(Exon *exon) {
  SliceAdaptor *sa;
  Slice *slice;
  AssemblyMapperAdaptor *ama;
  AssemblyMapper *assMapper;
  int sliceChrStart;
  int sliceChrEnd;
  int exonChrStart;
  int exonChrEnd;
  RawContigAdaptor *rca;
  MapperRangeSet *mapped;
  StringHash *SFHash = NULL;

  slice = (Slice *)Exon_getContig(exon);

  if (!slice) {
    fprintf(stderr,"Error: Cannot transform exon to raw contig unless it has an attached slice\n");
    exit(1);
  }

  sa = (SliceAdaptor *)Slice_getAdaptor(slice);

  if (!sa) {
    fprintf(stderr,"Error: Cannot transform exon to raw contig unless attached slice"
                   " has adaptor defined. (i.e. exon->contig->adaptor)\n");
    exit(1);
  }

  ama = DBAdaptor_getAssemblyMapperAdaptor(sa->dba);

  assMapper     = AssemblyMapperAdaptor_fetchByType(ama, CoordSystem_getName(Slice_getCoordSystem(slice)));
  rca           = DBAdaptor_getRawContigAdaptor(sa->dba);
  sliceChrStart = Slice_getChrStart(slice);
  sliceChrEnd   = Slice_getChrEnd(slice);


  if (Slice_getStrand(slice) == 1) {
    exonChrStart = Exon_getStart(exon) + sliceChrStart - 1;
    exonChrEnd   = Exon_getEnd(exon)   + sliceChrStart - 1;
  }
  else {
    exonChrEnd   = sliceChrEnd - Exon_getStart(exon) + 1;
    exonChrStart = sliceChrEnd - Exon_getEnd(exon)   + 1;
  }

  mapped = AssemblyMapper_mapCoordinatesToRawContig(assMapper,
     Slice_getChrName(slice),//Slice_getChrId(slice),
     exonChrStart,
     exonChrEnd,
     Exon_getStrand(exon)*Slice_getStrand(slice)
    );

  if (!mapped || !mapped->nRange) {
    fprintf(stderr, "Error: Exon couldnt map" );
    return exon;
  }

  // transform the supporting features to raw contig coords (hashed on contig)

// NOTE: Deliberately direct access 
  if (exon->supportingFeatures) {
    SFHash = StringHash_new(STRINGHASH_SMALL);
    Vector *sfs = Exon_getAllSupportingFeatures(exon);
    int i;
    for (i=0;i<Vector_getNumElement(sfs);i++) {
      BaseAlignFeature *baf = Vector_getElementAt(sfs,i);
      int j;
// NIY No eval
      Vector *mappedFeats = BaseAlignFeature_transformToRawContig(baf);
/*
      if ($@){
        $self->warn("Supporting feature didn't mapped ignoring $@");
        next SUPPORTING;
      }
*/
      for (j=0; j<Vector_getNumElement(mappedFeats); j++) {
        BaseAlignFeature *rcbaf = Vector_getElementAt(mappedFeats,j);
        char *contigName = BaseContig_getName(BaseAlignFeature_getSlice(rcbaf));
        Vector *rcVect;
        if (!StringHash_contains(SFHash, contigName)) {
          StringHash_add(SFHash, contigName, Vector_new());
        }
        rcVect = StringHash_getValue(SFHash,contigName);
        Vector_addElement(rcVect, rcbaf);
      }
    }
  }

    // thats a simple exon

    RawContig *rawContig;
    Exon *newExon;
    MapperRange *mr;
    MapperCoordinate *mc;
     
    mr = MapperRangeSet_getRangeAt(mapped, 0);

    if (mr->rangeType == MAPPERRANGE_GAP) {
      fprintf(stderr,"Error: exon lies on a gap cannot be mapped\n");
      exit(1);
    }

    mc = (MapperCoordinate *)mr;

    rawContig = RawContigAdaptor_fetchByDbID(rca, mc->id);
    newExon = Exon_new();

    // copy this exon
    Exon_copy(newExon, exon, SHALLOW_DEPTH);

    Exon_setStart(newExon, mc->start);
    Exon_setEnd(newExon, mc->end);
    Exon_setStrand(newExon, mc->strand);
    // attach raw contig
    Exon_setContig(newExon, rawContig);

    // replace old supporting feats with transformed supporting feats
    if (StringHash_contains(SFHash, RawContig_getName(rawContig))) {
      Vector *sfVect = StringHash_getValue(SFHash, RawContig_getName(rawContig));
      Exon_addSupportingFeatures(newExon, sfVect);
    }

    return newExon;
// NIY freeing old exons, old supporting features, SFHash Vectors
  if (SFHash) StringHash_free(SFHash, NULL);
}

void Exon_loadGenomicMapperImpl(Exon *exon, Mapper *mapper, IDType id, int start) {

// NIY Make the Exon_getContig consistent
  Mapper_addMapCoordinates( mapper, id, start, start+Exon_getLength(exon)-1,
                            Exon_getStrand(exon), BaseContig_getDbID(Exon_getContig(exon)),
                            Exon_getStart(exon),  Exon_getEnd(exon) );
}


void Exon_freeImpl(Exon *exon) {
  //printf("Exon_free not implemented\n");
  Object_decRefCount(exon);

  //fprintf(stderr,"Exon_free called refCount after dec = %d\n", Object_getRefCount(exon));
  if (Object_getRefCount(exon) > 0) {
    //fprintf(stderr,"Exon_free called BUT refCount after dec = %d so not freeing\n", Object_getRefCount(exon));
    return;
  } else if (Object_getRefCount(exon) < 0) {
    //fprintf(stderr,"Error: Negative reference count for Exon\n"
    //               "       Freeing it anyway\n");
  }

  //fprintf(stderr, "freeing exon %p %s nsupport = %d\n", exon, Exon_getStableId(exon), exon->supportingFeatures? Vector_getNumElement(exon->supportingFeatures) : 0);
  StableIdInfo_freePtrs(&(exon->si));
  if (Exon_getSeqCacheString(exon)!=NULL) {
    free(Exon_getSeqCacheString(exon));
  }

  if (exon->supportingFeatures) {
    Vector_setFreeFunc(exon->supportingFeatures, DNAAlignFeature_freeImpl);
    Vector_free(exon->supportingFeatures);
  }
  
  free(exon);
}

