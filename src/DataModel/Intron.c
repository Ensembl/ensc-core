/*
 * Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

/*
=head1 NAME Bio::EnsEMBL::Intron - A class representing an Intron
*/

#define __INTRON_MAIN__
#include "Intron.h"
#undef  __INTRON_MAIN__

#include "EnsC.h"

#include <string.h>

/*
=head2 new

  Arg [1]    : Bio::EnsEMBL::Exon The 5' exon for the intron; required
  Arg [2]    : Bio::EnsEMBL::Exon The 3' exon for the intron; required
  Arg [3]    : Bio::EnsEMBL::Analysis Analysis to link to this Intron
  Example    : $intron = new Bio::EnsEMBL::Intron($exon1, $exon2)
  Description: Create an Intron object from two exons and an optional analysis
  Returntype : Bio::EnsEMBL::Intron
  Exceptions : exons not on the same strand or slice.
  Caller     : general
  Status     : Stable

=cut
*/

Intron *Intron_new(Exon *e1, Exon *e2, Analysis *analysis) {
  Intron *intron;

  if ((intron = (Intron *)calloc(1,sizeof(Intron))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for intron\n");
    return NULL;
  }

  intron->objectType = CLASS_INTRON;
  Object_incRefCount(intron);

  intron->funcs = &intronFuncs;


  if (Exon_getStrand(e1) == -1 ) {
    Intron_setEnd(intron, Exon_getStart(e1) - 1);
    Intron_setStart(intron, Exon_getEnd(e2) + 1);
  } else {
    Intron_setStart(intron, Exon_getEnd(e1) + 1);
    Intron_setEnd(intron, Exon_getStart(e2) - 1);
  }

  if ( Exon_getStrand(e1) != Exon_getStrand(e2)) {
// Note perl had this commented out
    fprintf(stderr,"Exons on different strand. Not allowed\n");
    exit(1);
  } else {
    Intron_setStrand(intron, Exon_getStrand(e1));
  }

  if (Exon_getSlice(e1) != Exon_getSlice(e2)) {
// Perl had a wierd condition with &&s rather than ||s. Didn't look right
    if ( ( strcmp(Slice_getSeqRegionName(Exon_getSlice(e1)), Slice_getSeqRegionName(Exon_getSlice(e2)))) || 
         ( strcmp(Slice_getSeqRegionName(Exon_getSlice(e1)), Slice_getSeqRegionName(Exon_getSlice(e2))))) {
      fprintf(stderr,"Intron with Exons on different slices. Not allowed\n");
      exit(1);
    } else {
      fprintf(stderr, "Warning: Exons have different slice references to the same seq_region\n");
// Do we still want to set slice???
    }
  } else {
    Intron_setSlice(intron, Exon_getSlice(e1));
  }
  
  if (analysis != NULL) {
    Intron_setAnalysis(intron, analysis);
  }

  Intron_setPrevExon(intron, e1);
  Intron_setNextExon(intron, e2);

  return intron;
}

/*
=head2 length

  Args       : none
  Example    : $length = $intron->length();
  Description: Returns the length of this intron
  Returntype : Integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/
long Intron_getLength(Intron *intron) {

  return Intron_getEnd(intron) - Intron_getStart(intron) + 1;
}

/*
=head2 is_splice_canonical

  Example     : my $canonical = $intron->is_splice_canonical(); 
  Description : Indicates if the splice site is considered normal. This means
                splice site variants equal to (D == donor, A == acceptor)
                  GT (D) => AG (A) 
                  AT (D) => AC (A)
                  GC (D) => AG (A)
  Returntype  : Boolean indicating if the splice was as expected
  Exceptions  : See splice_seq

=cut
*/
static int warnedAboutCanon = 0;
int Intron_getIsSpliceCanonical(Intron *intron) {

  if (!warnedAboutCanon) {
    fprintf(stderr,"NIY: Intron_getIsSpliceCanonical\n");
    warnedAboutCanon = 1;
  }
/* NIY
  my $splice = join q{}, @{$self->splice_seq()}; 
  my $canonical = {
    'GTAG' => 1, 'ATAC' => 1, 'GCAG' => 1
  }->{$splice};
  return $canonical || 0;
*/
  return 1;
}

/*
=head2 splice_seq

  Example     : my ($donor, $acceptor) = @{$intron->splice_seq}; 
  Description : Get the donor and acceptor splice sites for this intron
  Returntype  : ArrayRef[String] The donor and acceptor sequences as Strings
  Exceptions  : Thrown if a feature Slice cannot be found

=cut
*/
/* NIY
sub splice_seq {
  my ($self) = @_;
  my $slice = $self->feature_Slice();
  throw "Cannot retrieve feature_Slice() for this Intron" unless $slice;
  my $length = $self->length();
  my $donor_seq    = uc($slice->subseq(1,2));
  my $acceptor_seq = uc($slice->subseq($length - 1, $length));
  return [$donor_seq, $acceptor_seq];
}
*/

void Intron_freeImpl(Intron *intron) {
  //fprintf(stderr,"NIY: Intron_freeImpl\n");
  free(intron);
}
