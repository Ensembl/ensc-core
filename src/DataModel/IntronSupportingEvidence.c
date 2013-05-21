/*
=head1 DESCRIPTION

Formalises an Intron with information about why it is a believed Intron. This 
serves as a parallel object to Bio::EnsEMBL::Intron which you can use 
to populate values in this field from. They are different objects though
due to Intron's non-existence as a DB data structure.
*/

#define __INTRONSUPPORTINGEVIDENCE_MAIN__
#include "IntronSupportingEvidence.h"
#undef  __INTRONSUPPORTINGEVIDENCE_MAIN__
#include "IntronSupportingEvidenceAdaptor.h"
#include "Intron.h"
#include "StrUtil.h"
#include "Transcript.h"

#include "EnsC.h"

#include <strings.h>


static char *supportedTypes[] = { "NONE" , "DEPTH", NULL };

/*
=head2 new

  Arg [-ANALYSIS]     : Bio::EnsEMBL::Analysis The analysis this intron is linked to
  Arg [-START]        : int - start postion of the IntronSupportingEvidence
  Arg [-END]          : int - end position of the IntronSupportingEvidence
  Arg [-STRAND]       : int - strand the IntronSupportingEvidence is on
  Arg [-SLICE]        : Bio::EnsEMBL::Slice - the slice the IntronSupportingEvidence is on
  Arg [-INTRON]       : Bio::EnsEMBL::Intron Intron the evidence is based 
                        on. Useful if you are not specifying the location 
                        parameters as we will take them from this 
  Arg [-HIT_NAME]     : String The name of the hit
  Arg [-SCORE]        : Double The score associated with the supporting evidence
  Arg [-SCORE_TYPE]   : String The type of score we are representing
  Example             : Bio::EnsEMBL::IntronSupportingEvidence->new(
                          -ANALYSIS => $analysis, -INTRON => $intron, 
                          -SCORE => 100, -SCORE_TYPE => 'DEPTH');
  Description         : Returns a new instance of this object
  Returntype          : Bio::EnsEMBL::IntronSupportEvidence
  Exceptions          : Thrown if data is not as requested

=cut
*/
IntronSupportingEvidence *IntronSupportingEvidence_new(void) {
  IntronSupportingEvidence *ise;

  if ((ise = (IntronSupportingEvidence *)calloc(1,sizeof(IntronSupportingEvidence))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for IntronSupportingEvidence\n");
    return NULL;
  }

  ise->objectType = CLASS_INTRONSUPPORTINGEVIDENCE;
  Object_incRefCount(ise);

  ise->funcs = &intronSupportingEvidenceFuncs;

  return ise;
}


/* Note does some stuff if given an intron
sub new {
  my ($class, @args) = @_;
  
  my $self = $class->SUPER::new(@args);
  
  my ($intron, $hit_name, $score, $score_type, $is_splice_canonical) = 
    rearrange([qw/intron hit_name score score_type is_splice_canonical/], @args);
  
  if($intron) {
    $self->set_values_from_Intron($intron);
  }
  $self->hit_name($hit_name) if $hit_name;
  $self->score($score) if $score;
  $self->score_type($score_type) if $score_type;
  $self->is_splice_canonical($is_splice_canonical) if $is_splice_canonical;
  
  return $self;
}
*/

/*
=head2 set_values_from_Intron

  Arg [1]     : Bio::EnsEMBL::Intron The intron to base this object on
  Example     : $ise->set_values_from_Intron($intron);
  Description : Sets the start, end, strand and slice of this ISE instance
                using values from the given Intron object.
  Returntype  : None
  Exceptions  : Thrown if data is not as requested

=cut
*/
void IntronSupportingEvidence_setValuesFromIntron(IntronSupportingEvidence *ise, Intron *intron) {

  IntronSupportingEvidence_setStart (ise, Intron_getStart(intron));
  IntronSupportingEvidence_setEnd   (ise, Intron_getEnd(intron));
  IntronSupportingEvidence_setStrand(ise, Intron_getStrand(intron));
  IntronSupportingEvidence_setSlice (ise, Intron_getSlice(intron));
  IntronSupportingEvidence_setIsSpliceCanonical(ise, Intron_getIsSpliceCanonical(intron));

  return;
}

/*
=head2 get_Intron

  Arg [1]     : Bio::EnsEMBL::Transcript
  Example     : my $intron = $ise->intron($transcript);
  Description : Provides access to an Intron object by using a given transcript 
                object and its associcated array of Exons.
  Returntype  : Bio::EnsEMBL::Intron
  Exceptions  : None

=cut
*/
Intron *IntronSupportingEvidence_getIntron(IntronSupportingEvidence *ise, Transcript *transcript) {

  Exon *fivePrime  = IntronSupportingEvidence_findPreviousExon(ise, transcript);
  Exon *threePrime = IntronSupportingEvidence_findNextExon(ise, transcript);
  
  Intron *intron = Intron_new(fivePrime, threePrime, NULL);
  
  return intron;
}

/*
=head2 hit_name

  Arg [1]     : String name of the hit
  Example     : $ise->hit_name('hit');
  Description : Getter/setter for hit name i.e. an identifier for the alignments
  Returntype  : String
  Exceptions  : None

=cut
*/
char *IntronSupportingEvidence_setHitName(IntronSupportingEvidence *ise, char *hitName) {
  StrUtil_copyString(&(ise->hitName), hitName, 0);

  return ise->hitName;
}

/*
=head2 score_type

  Arg [1]     : String the enum type. Currently only allowed NONE or DEPTH
  Example     : $ise->score_type('DEPTH');
  Description : Gets and sets the type of score this instance represents
  Returntype  : String
  Exceptions  : Thrown if given an unsupported type of data

=cut
*/
char *IntronSupportingEvidence_setScoreType(IntronSupportingEvidence *ise, char *scoreType) {
  if (scoreType == NULL) {
    fprintf(stderr,"NULL score type in ISE_setScoreType\n");
    exit(1);
  }

  int found = 0;
  int i;
  for (i=0; supportedTypes[i]!=NULL && !found; i++) {
    if (!strcasecmp(supportedTypes[i], scoreType)) {
      found = 1;
    }
  }
  if (!found) {
    fprintf(stderr, "The score type '%s' is not allowed. Allowed values are [", scoreType);
    for (i=0; supportedTypes[i]!=NULL; i++) {
      if (i!= 0) {
        fprintf(stderr,", %s",supportedTypes[i]); 
      } else {
        fprintf(stderr,"%s",supportedTypes[i]); 
      }
    }
    fprintf(stderr,"]\n");
    exit(1);
  }
   
  EcoString_copyStr(ecoSTable, &(ise->scoreType), scoreType, 0);

  if (ise->scoreType == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for scoreType\n");
    return NULL;
  }

  return ise->scoreType;
}

/*
=head2 has_linked_transcripts

  Example     : $ise->has_linked_transcripts();
  Description : Returns true if we have transcripts linked to this ISE
  Returntype  : Boolean
  Exceptions  : Thrown if we do not have an attached adaptor

=cut
*/
int IntronSupportingEvidence_hasLinkedTranscripts(IntronSupportingEvidence *ise) {
  if (IntronSupportingEvidence_getAdaptor(ise) == NULL) {
    fprintf(stderr,"No attached adaptor. Cannot find linked Transcripts unless this is a persistent object\n");
    exit(1);
  }

  IntronSupportingEvidenceAdaptor *isea = (IntronSupportingEvidenceAdaptor *)IntronSupportingEvidence_getAdaptor(ise);

  Vector *transcriptIds = IntronSupportingEvidenceAdaptor_listLinkedTranscriptIds(isea, ise);

  int hasLinked = (Vector_getNumElement(transcriptIds) > 0);

  Vector_free(transcriptIds);

  return hasLinked;
}

/* NIY - maybe never as I can't see when this is used
=head2 equals

  Arg [1]     : Bio::EnsEMBL::IntronSupportEvidence Object to compare to
  Example     : $ise->equals($another_ise);
  Description : Asserts if the given IntronSupportEvidence instance was equal to this
  Returntype  : Boolean
  Exceptions  : None

=cut

sub equals {
  my ($self, $other) = @_;
  my $equal = $self->SUPER::equals($other);
  return 0 if ! $equal;
  return ( 
    ($self->hit_name()||q{}) eq ($other->hit_name()||q{}) &&
    ($self->score_type() eq $other->score_type()) &&
    ($self->score() == $other->score())) ? 1 : 0;
}
*/

/*
=head2 find_previous_Exon

  Arg [1]     : Bio::EnsEMBL::Transcript Transcript to search for the Exons from
  Example     : $ise->find_previous_Exon($transcript);
  Description : Loops through those Exons available from the Transcript and
                attempts to find one which was the 5' flanking exon. If the
                object has already been persisted we will use dbIDs to
                find the Exons
  Returntype  : Bio::EnsEMBL::Exon
  Exceptions  : None

=cut
*/
Exon *IntronSupportingEvidence_findPreviousExon(IntronSupportingEvidence *ise, Transcript *transcript) {
  
  //Use DB IDs if we have them
  IDType exonId = 0;

  if (IntronSupportingEvidence_getAdaptor(ise) == NULL) {

    IntronSupportingEvidenceAdaptor *isea = (IntronSupportingEvidenceAdaptor *)IntronSupportingEvidence_getAdaptor(ise);

    IDType flankingIds[2];
    
    if (IntronSupportingEvidenceAdaptor_fetchFlankingExonIds(isea, ise, transcript, flankingIds) != NULL) {
      exonId = flankingIds[0];
    
    } else {
      fprintf(stderr,"Odd - findPreviousExon with adaptor but didn't find prev exon id in db\n");
    }
  }
  
  long start = IntronSupportingEvidence_getStart(ise);
  long end   = IntronSupportingEvidence_getEnd(ise);
  int strand = IntronSupportingEvidence_getStrand(ise);

  int i;
  for (i=0; i<Transcript_getExonCount(transcript); i++) {
    Exon *exon = Transcript_getExonAt(transcript, i);
    if (exonId) {
      if (exonId == Exon_getDbID(exon)) {
        return exon;
      }
    } else {
      if(strand == 1) {
        if (Exon_getEnd(exon) == start-1) {
          return exon;
        }
      } else {
        if (Exon_getStart(exon) == end+1) {
          return exon;
        }
      }
    }
  }
  return NULL;
}

/*
=head2 find_next_Exon

  Arg [1]     : Bio::EnsEMBL::Transcript Transcript to search for the Exons from
  Example     : $ise->find_next_Exon($transcript);
  Description : Loops through those Exons available from the Transcript and
                attempts to find one which was the 3' flanking exon. If the
                object has already been persisted we will use dbIDs to
                find the Exons
  Returntype  : Bio::EnsEMBL::Exon
  Exceptions  : None

=cut
*/
Exon *IntronSupportingEvidence_findNextExon(IntronSupportingEvidence *ise, Transcript *transcript) {
  
  //Use DB IDs if we have them
  IDType exonId = 0;

  if (IntronSupportingEvidence_getAdaptor(ise) == NULL) {

    IntronSupportingEvidenceAdaptor *isea = (IntronSupportingEvidenceAdaptor *)IntronSupportingEvidence_getAdaptor(ise);

    IDType flankingIds[2];
    
    if (IntronSupportingEvidenceAdaptor_fetchFlankingExonIds(isea, ise, transcript, flankingIds) != NULL) {
      exonId = flankingIds[1];
    
    } else {
      fprintf(stderr,"Odd - findNextExon with adaptor but didn't find prev exon id in db\n");
    }
  }
  
  long start = IntronSupportingEvidence_getStart(ise);
  long end   = IntronSupportingEvidence_getEnd(ise);
  int strand = IntronSupportingEvidence_getStrand(ise);

  int i;
  for (i=0; i<Transcript_getExonCount(transcript); i++) {
    Exon *exon = Transcript_getExonAt(transcript, i);
    if (exonId) {
      if (exonId == Exon_getDbID(exon)) {
        return exon;
      }
    } else {
      if(strand == 1) {
        if (Exon_getStart(exon) == end+1) {
          return exon;
        }
      } else {
        if (Exon_getEnd(exon) == start-1) {
          return exon;
        }
      }
    }
  }
  return NULL;
}


void IntronSupportingEvidence_freeImpl(IntronSupportingEvidence *ise) {
  fprintf(stderr, "NIY: IntronSupportingEvidenceFree\n");

}
