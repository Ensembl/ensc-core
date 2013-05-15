#define __SEQFEATURE_MAIN__
#include "SeqFeature.h"
#undef __SEQFEATURE_MAIN__

#include "DBAdaptor.h"
#include "StrUtil.h"
#include "AssemblyMapperAdaptor.h"
#include "RawContigAdaptor.h"
#include "SliceAdaptor.h"
#include "SeqFeatureFactory.h"
#include "CoordSystemAdaptor.h"
#include "SimpleFeature.h"
#include "ProjectionSegment.h"

SeqFeature *SeqFeature_new(void) {
  SeqFeature *sf;

  if ((sf = (SeqFeature *)calloc(1,sizeof(SeqFeature))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for seq feature\n");
    return NULL;
  }

  sf->funcs = &seqFeatureFuncs;

  return sf;
}


ECOSTRING SeqFeature_setSeqName(SeqFeature *sf, char *seqName) {
  EcoString_copyStr(ecoSTable, &(sf->seqName),seqName,0);

  return sf->seqName;
}

ECOSTRING SeqFeature_getSeqName(SeqFeature *sf) {
  BaseContig *contig = SeqFeature_getContig(sf);

  if (contig) {
    return BaseContig_getName(contig);
  } else if (sf->seqName) {
    return sf->seqName;
  } else {
    fprintf(stderr,"Warning: No seq name defined for feature\n");
    return emptyString;
  }
}

int SeqFeature_startCompFunc(const void *a, const void *b) {
  SeqFeature **e1 = (SeqFeature **)a;
  SeqFeature **e2 = (SeqFeature **)b;

  if (SeqFeature_getStart(*e1) > SeqFeature_getStart(*e2)) {
    return 1;
  } else if (SeqFeature_getStart(*e1) < SeqFeature_getStart(*e2)) {
    return -1;
  } else {
    return 0;
  }
}

int SeqFeature_reverseStartCompFunc(const void *a, const void *b) {
  SeqFeature **e1 = (SeqFeature **)a;
  SeqFeature **e2 = (SeqFeature **)b;

  if (SeqFeature_getStart(*e2) > SeqFeature_getStart(*e1)) {
    return 1;
  } else if (SeqFeature_getStart(*e2) < SeqFeature_getStart(*e1)) {
    return -1;
  } else {
    return 0;
  }
}

void SeqFeature_freePtrs(SeqFeature *sf) {
  if (sf->seqName)  EcoString_freeStr(ecoSTable, sf->seqName);
  if (sf->analysis) Analysis_free(sf->analysis);
// NIY Is this the right thing to do
  // NIY Freeing contig if (sf->contig)   BaseContig_free(sf->contig);
}

/* Note this code comes from Feature.pm not SeqFeature.pm
   For now I'm only interested in implementing the transform and transfer methods
   I'll get to the others later

=head2 new

  Arg [-SLICE]: Bio::EnsEMBL::SLice - Represents the sequence that this
                feature is on. The coordinates of the created feature are
                relative to the start of the slice.
  Arg [-START]: The start coordinate of this feature relative to the start
                of the slice it is sitting on.  Coordinates start at 1 and
                are inclusive.
  Arg [-END]  : The end coordinate of this feature relative to the start of
                the slice it is sitting on.  Coordinates start at 1 and are
                inclusive.
  Arg [-STRAND]: The orientation of this feature.  Valid values are 1,-1,0.
  Arg [-SEQNAME] : A seqname to be used instead of the default name of the 
                of the slice.  Useful for features that do not have an 
                attached slice such as protein features.
  Arg [-dbID]   : (optional) internal database id
  Arg [-ADAPTOR]: (optional) Bio::EnsEMBL::DBSQL::BaseAdaptor
  Example    : $feature = Bio::EnsEMBL::Feature->new(-start    => 1, 
                                                     -end      => 100,
                                                     -strand   => 1,
                                                     -slice    => $slice,
                                                     -analysis => $analysis);
  Description: Constructs a new Bio::EnsEMBL::Feature.  Generally subclasses
               of this method are instantiated, rather than this class itself.
  Returntype : Bio::EnsEMBL::Feature
  Exceptions : Thrown on invalid -SLICE, -ANALYSIS, -STRAND ,-ADAPTOR arguments
  Caller     : general, subclass constructors
  Status     : Stable

=cut


sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;
  my ( $start, $end, $strand, $slice, $analysis,$seqname, $dbID, $adaptor ) =
      rearrange(['START','END','STRAND','SLICE','ANALYSIS', 'SEQNAME',
		 'DBID', 'ADAPTOR'], @_);   
  if($slice) {
    if(!ref($slice) || !($slice->isa('Bio::EnsEMBL::Slice') or $slice->isa('Bio::EnsEMBL::LRGSlice')) ) {
      throw('-SLICE argument must be a Bio::EnsEMBL::Slice not '.$slice);
    }
  }

  if($analysis) {
    if(!ref($analysis) || !$analysis->isa('Bio::EnsEMBL::Analysis')) {
      throw('-ANALYSIS argument must be a Bio::EnsEMBL::Analysis not '.
            $analysis);
    }
  }

  if(defined($strand)) {
    if(!($strand == 1) && !($strand == -1) && !($strand == 0)) {
      throw('-STRAND argument must be 1, -1, or 0');
    }
  }

  if(defined($start) && defined($end)) {
      if (($start =~ /\d+/) && ($end =~ /\d+/)) {
	  if($end+1 < $start) {
	      throw(sprintf('Start (%d) must be less than or equal to end+1 (%d)', $start, ($end+1)));
	  }
      } else {
	      throw('Start and end must be integers');
      }
  }

  my $self =  bless({'start'    => $start,
                'end'      => $end,
                'strand'   => $strand,
                'slice'    => $slice,
                'analysis' => $analysis,
                'seqname'  => $seqname,
                'dbID'     => $dbID}, $class);

  $self->adaptor($adaptor);
  return $self;
}


=head2 new_fast

  Arg [1]    : hashref to be blessed
  Description: Construct a new Bio::EnsEMBL::Feature using the hashref.
  Exceptions : none
  Returntype : Bio::EnsEMBL::Feature
  Caller     : general, subclass constructors
  Status     : Stable

=cut


sub new_fast {
  my $class = shift;
  my $hashref = shift;
  my $self = bless $hashref, $class;
  weaken($self->{adaptor})  if ( ! isweak($self->{adaptor}) );
  return $self;
}

=head2 start

  Arg [1]    : (optional) int $start
               The start of this feature relative to the start of the slice
               that it is on.
  Example    : $start = $feat->start()
  Description: Getter/Setter for the start of this feature relative to the 
               start of the slice it is on.  Note that negative values, or
               values exceeding the length of the slice are permitted.
               Start must be less than or equal to the end regardless of the 
               strand. Coordinate values start at 1 and are inclusive.
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub start {
  my ( $self, $value ) = @_;

  if ( defined($value) ) {
    $self->{'start'} = $value;
  }

  return $self->{'start'};
}



=head2 end

  Arg [1]    : (optional) int $end
  Example    : $end = $feat->end();
  Description: Getter/Setter for the end of this feature relative to the
               start of the slice that it is on.  Note that negative values,
               of values exceeding the length of the slice are permitted.  End
               must be greater than or equal to start regardless of the strand.
               Coordinate values start at 1 and are inclusive.
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub end {
  my ( $self, $value ) = @_;

  if ( defined($value) ) {
    $self->{'end'} = $value;
  }

  return $self->{'end'};
}




=head2 strand

  Arg [1]    : (optional) int $strand
  Example    : $feat->strand(-1);
  Description: Getter/Setter for the strand of this feature relative to the
               slice it is on.  0 is an unknown or non-applicable strand.  
               -1 is the reverse (negative) strand and 1 is the forward 
               (positive) strand.  No other values are permitted.
  Returntype : int
  Exceptions : thrown if an invalid strand argument is passed
  Caller     : general
  Status     : Stable

=cut

sub strand {
  my ( $self, $strand ) = @_;

  if ( defined($strand) ) {
    if ( $strand != 0 && $strand != 1 && $strand != -1 ) {
      throw('strand argument must be 0, -1 or 1');
    }

    $self->{'strand'} = $strand;
  }

  return $self->{'strand'};
}

=head2 move

  Arg [1]    : int start
  Arg [2]    : int end
  Arg [3]    : (optional) int strand
  Description: Sets the start, end and strand in one call rather than in 
               3 seperate calls to the start(), end() and strand() methods.
               This is for convenience and for speed when this needs to be
               done within a tight loop.
  Returntype : none
  Exceptions : Thrown is invalid arguments are provided
  Caller     : general
  Status     : Stable

=cut

sub move {
  my $self = shift;

  throw('start and end arguments are required') if(@_ < 2);

  my $start  = shift;
  my $end    = shift;
  my $strand = shift;

  if(defined($start) && defined($end) && $end < $start) {
    throw('start must be less than or equal to end');
  }
  if(defined($strand) && $strand != 0 && $strand != -1 && $strand != 1) {
    throw('strand must be 0, -1 or 1');
  }

  $self->{'start'} = $start;
  $self->{'end'} = $end;
  $self->{'strand'} = $strand if(defined($strand));
}



=head2 length

  Arg [1]    : none
  Example    : $length = $feat->length();
  Description: Returns the length of this feature
  Returntype : Integer
  Exceptions : Throws if end < start and the feature is not on a
               circular slice
  Caller     : general
  Status     : Stable

=cut

sub length {
  my ($self) = @_;

  if ( $self->{'end'} < $self->{'start'} ) {
    # if circular, we can work out the length of an origin-spanning
    # feature using the size of the underlying region.
    if ( $self->slice() && $self->slice()->is_circular() ) {
      my $len =
        $self->slice()->seq_region_length() -
        ( $self->{'start'} - $self->{'end'} ) + 1;
      return $len;
    } else {
      throw(   "Cannot determine length of non-circular feature "
             . "where start > end" );
    }
  }

  return $self->{'end'} - $self->{'start'} + 1;
}

=head2 analysis

  Arg [1]    : (optional) Bio::EnsEMBL::Analysis $analysis
  Example    : $feature->analysis(new Bio::EnsEMBL::Analysis(...))
  Description: Getter/Setter for the analysis that is associated with 
               this feature.  The analysis describes how this feature 
               was derived.
  Returntype : Bio::EnsEMBL::Analysis
  Exceptions : thrown if an invalid argument is passed
  Caller     : general
  Status     : Stable

=cut

sub analysis {
  my $self = shift;

  if(@_) {
    my $an = shift;
    if(defined($an) && (!ref($an) || !$an->isa('Bio::EnsEMBL::Analysis'))) {
      throw('analysis argument must be a Bio::EnsEMBL::Analysis');
    }
    $self->{'analysis'} = $an;
  }

  return $self->{'analysis'};
}



=head2 slice

  Arg [1]    : (optional) Bio::EnsEMBL::Slice $slice
  Example    : $seqname = $feature->slice()->name();
  Description: Getter/Setter for the Slice that is associated with this 
               feature.  The slice represents the underlying sequence that this
               feature is on.  Note that this method call is analagous to the
               old SeqFeature methods contig(), entire_seq(), attach_seq(),
               etc.
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : thrown if an invalid argument is passed
  Caller     : general
  Status     : Stable

=cut

sub slice {
  my ( $self, $slice ) = @_;

  if ( defined($slice) ) {
    if (    !check_ref( $slice, 'Bio::EnsEMBL::Slice' )
         && !check_ref( $slice, 'Bio::EnsEMBL::LRGSlice' ) )
    {
      throw('slice argument must be a Bio::EnsEMBL::Slice');
    }

    $self->{'slice'} = $slice;
  } elsif ( @_ > 1 ) {
    delete($self->{'slice'});
  }

  return $self->{'slice'};
}

=head2 equals

  Arg [1]       : Bio::EnsEMBL::Feature object
  Example       : if ($featureA->equals($featureB)) { ... }
  Description   : Compares two features using various criteria.  The
                  test for eqality goes through the following list and
                  terminates at the first true match:

                  1. If the two features are the same object, they are
                     equal.
                  2. If they are of different types (e.g., transcript
                     and gene), they are *not* equal.
                  3. If they both have dbIDs: if these are the same,
                     then they are equal, otherwise not.
                  4. If they both have slices and analysis objects:
                     if the analysis dbIDs are the same and the
                     seq_region_id are the same, along with
                     seq_region_start and seq_region_end, then they are
                     equal, otherwise not.

                  If none of the above is able to determine equality,
                  undef is returned.

    Return type : tri-Boolean (0, 1, undef = "unknown")

    Exceptions  : Thrown if a non-feature is passed as the argument.

=cut

sub equals {
  my ( $self, $feature ) = @_;

  # If the features are the same object, they are equal.
  if ( !defined($feature) ) { return 0 }
  if ( $self eq $feature ) { return 1 }

  assert_ref( $feature, 'Bio::EnsEMBL::Feature' );

  # If the features have different types, they are *not* equal.
  if ( ref($self) ne ref($feature) ) {
    return 0;
  }

  # If the features has the same dbID, they are equal.
  if ( defined( $self->dbID() ) && defined( $feature->dbID() ) ) {
    if   ( $self->dbID() == $feature->dbID() ) { return 1 }
    else                                       { return 0 }
  }

  # We now know that one of the features do not have a dbID.

  # If the features have the same start, end, strand and seq_region_id,
  # and analysis_id, they are equal.
  if (
     ( defined( $self->analysis() ) && defined( $feature->analysis() ) )
     && ( defined( $self->slice() ) && defined( $feature->slice() ) ) )
  {
    if ( ( $self->start() == $feature->start() ) &&
         ( $self->end() == $feature->end() ) &&
         ( $self->strand() == $feature->strand() ) &&
         ( $self->slice()->get_seq_region_id() ==
           $feature->slice()->get_seq_region_id() ) &&
         ( $self->analysis()->dbID() == $feature->analysis()->dbID() ) )
    {
      return 1;
    }
    else { return 0 }
  }

  # We now know that one of the features does not have either analysis
  # or slice.

  # We don't know if the features are equal.  This happens if they are
  # not the same object but are of the same type, and one of them lacks
  # dbID, and if there aren't slice and analysis objects attached to
  # them both.
  return undef;
} ## end sub equals
*/


/*
=head2 transform

  Arg [1]    : string $coord_system
               The coord system to transform this feature to.
  Arg [2]    : string $version (optional)
               The version of the coord system to transform this feature to.
  Arg [3]    : Bio::EnsEMBL::Slice (optional)
               Specified when a projection may land on many overlapping slices
               and disambiguation is required.
  Example    : $feature = $feature->transform('contig');
               next if(!defined($feature));
  Description: Returns a copy of this feature, but converted to a different
               coordinate system. The converted feature will be placed on a
               slice which spans an entire sequence region of the new
               coordinate system. If the requested coordinate system is the
               same coordinate system it is simply placed on a slice which
               spans the entire seq_region (as opposed to the original slice
               which may have only partially covered the seq_region).

               If a feature spans a boundary in the new coordinate system,
               undef is returned instead.

               For example, transforming an exon in contig coordinates to one 
               in chromosomal coodinates will place the exon on a slice of an 
               entire chromosome.
  Returntype : Bio::EnsEMBL::Feature (or undef)
  Exceptions : thrown if an invalid coordinate system is provided
               warning if Feature is not attached to a slice
  Caller     : general, transfer()
  Status     : Stable

=cut
*/
SeqFeature *SeqFeature_transform(SeqFeature *sf, char *csName, char *csVersion, Slice *toSlice) {
  // 
  // For backwards compatibility check if the arguments are old style args
  // 
  // Can't really do ref(cs_name) in C, so won't check for it
  if (csName == NULL) { //|| ref($cs_name)) {
    fprintf(stderr, "Deprecated: Calling transform without a coord system name is deprecated. Steve hasn't implemented the deprecated method so bye!");
    exit(1);
// Perl has    return $self->_deprecated_transform($cs_name);
  }

  Slice *slice = SeqFeature_getSlice(sf);

  if (slice == NULL) {
    fprintf(stderr, "warning: Feature cannot be transformed without attached slice.\n");
    return NULL;
  }

  if (!Slice_getAdaptor(slice)) {
    fprintf(stderr, "Feature cannot be transformed without adaptor on attached slice.\n");
    return NULL;
  }

  //use db from slice since this feature may not yet be stored in a database
  SliceAdaptor *sa = (SliceAdaptor *)Slice_getAdaptor(slice);
  CoordSystemAdaptor *csa = DBAdaptor_getCoordSystemAdaptor(sa->dba);

  CoordSystem *cs = CoordSystemAdaptor_fetchByName(csa, csName, csVersion);

  CoordSystem *currentCs = Slice_getCoordSystem(slice);

  if (currentCs == NULL) {
    fprintf(stderr, "Warning: Feature cannot be transformed without CoordSystem on attached slice.\n");
    return NULL;
  }

  if (cs == NULL) {
    fprintf(stderr, "Cannot transform to unknown coordinate system [%s %s]\n", csName, csVersion);
    exit(1);
  }

  // if feature is already in the requested coordinate system, we can just
  // return a copy
  if (!CoordSystem_compare(cs, currentCs) && 
      Slice_getStart(slice) == 1 &&
      Slice_getStrand(slice) == 1) {
    // Temporary Hack - really need a clone method
    SeqFeature *newFeature = SeqFeatureFactory_newFeatureFromFeature(sf);
    SeqFeature_setStart(newFeature, SeqFeature_getStart(sf));
    SeqFeature_setEnd(newFeature, SeqFeature_getEnd(sf));
    SeqFeature_setStrand(newFeature, SeqFeature_getStrand(sf));
    SeqFeature_setSlice(newFeature, SeqFeature_getSlice(sf));

    SeqFeature_setDbID(newFeature, SeqFeature_getDbID(sf));

    return newFeature;

//    my $new_feature;
//    %$new_feature = %$self;
//    bless $new_feature, ref $self;
//    return $new_feature;
  }

  Vector *projection;
  if (toSlice != NULL) {
    projection = SeqFeature_projectToSlice(sf, toSlice);
  } else{
    projection = SeqFeature_project(sf, csName, csVersion);
  }

  int nProjection = Vector_getNumElement(projection);

  if (nProjection == 0) {
    return NULL;
  }

  if (nProjection != 1 && toSlice == NULL) {
// Warns were commented out - I've reinstated them for now for C just in case they catch something
    //fprintf(stderr, "MORE than one projection and NO slice specified from %s to %s, %s\n", 
    //        Slice_getName(SeqFeature_getSlice(sf)), csName, csVersion);
    return NULL;
  }

  int index = 0;
  if (toSlice != NULL) {
    int found = 0;
    int i;

    for (i=0; i<nProjection; i++) {
      ProjectionSegment *proj = Vector_getElementAt(projection, i);
      Slice *projSlice = ProjectionSegment_getToSlice(proj);
// NIY: Here we use ids, in transfer we use names - which should it be????
//      Note was also using eq for ids - bit odd
      if (Slice_getSeqRegionId(toSlice) == Slice_getSeqRegionId(projSlice)){
	found = 1;
	index = i;
      }
    }
    if (! found) {
      if (nProjection != 1) {
	if (nProjection == 0) {
/* Code can never be reached because of check for 0 number of projections above!
	  warn "number of mappings is ".@$projection."\n";
	  warn "could not project feature ".ref($self)." from ".$self->slice->seq_region_name." to ".$to_slice->seq_region_name."\n";
	  warn "In the region of ".$self->slice->start." <-> ".$self->slice->end."\n";
	  warn "feat start=".($self->slice->start+$self->start)."\tend=".($self->slice->start+$self->end)."\n";
*/
	} else {
          for (i=0; i<nProjection; i++) {
            ProjectionSegment *proj = Vector_getElementAt(projection, i);
	    Slice *projSlice = ProjectionSegment_getToSlice(proj);
	    fprintf(stderr, "available slice %s\n", Slice_getSeqRegionName(projSlice));
	  }
	  fprintf(stderr, "Warning: MORE than one projection and none to slice specified (%s)\n", Slice_getSeqRegionName(toSlice));
	}
      }	else {
        for (i=0; i<nProjection; i++) {
          ProjectionSegment *proj = Vector_getElementAt(projection, i);
	  Slice *projSlice = ProjectionSegment_getToSlice(proj);
	  fprintf(stderr, "Mapping is to %s\n", Slice_getSeqRegionName(projSlice));
	}
	fprintf(stderr,"One projection but none to slice specified\n");
      }
      return NULL;
    }
  }
 
  Slice *pSlice = ProjectionSegment_getToSlice((ProjectionSegment *)Vector_getElementAt(projection, index));

  CoordSystem *pCs = Slice_getCoordSystem(pSlice);

  Slice *newSlice = SliceAdaptor_fetchByRegion(sa, 
                                               CoordSystem_getName(pCs),
					       Slice_getSeqRegionName(pSlice),
					       POS_UNDEF, // start
					       POS_UNDEF, // end
					       1,         // strand
					       CoordSystem_getVersion(pCs),
                                               0);
  
  SeqFeature *newFeature;
//NIY: How to do the copy - that is the question??
//  %$new_feature = %$self;
//  bless $new_feature, ref $self;
// Temporary hack
  newFeature = SeqFeatureFactory_newFeatureFromFeature(sf);

  SeqFeature_setStart(newFeature, Slice_getStart(pSlice));
  SeqFeature_setEnd(newFeature, Slice_getEnd(pSlice));
      
  SeqFeature_setDbID(newFeature, SeqFeature_getDbID(sf));

// Huhhh?  What's this strand 0 stuff???????
  if (SeqFeature_getSlice(sf) == 0) {
    SeqFeature_setStrand(newFeature, 0);
  } else {
    SeqFeature_setStrand(newFeature, Slice_getStrand(pSlice));
  }

  SeqFeature_setSlice(newFeature, newSlice);

  // NIY: Free stuff (projection segments etc)

  return newFeature;
}

/*
=head2 transfer

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               The slice to transfer this feature to
  Example    : $feature = $feature->transfer($slice);
               next if(!defined($feature));
  Description: Returns a copy of this feature which has been shifted onto
               another slice.

               If the new slice is in a different coordinate system the
               feature is transformed first and then placed on the slice.
               If the feature would be split across a coordinate system
               boundary or mapped to a gap undef is returned instead.

               If the feature cannot be placed on the provided slice because
               it maps to an entirely different location, undef is returned
               instead.

  Returntype : Bio::EnsEMBL::Feature (or undef)
  Exceptions : throw on incorrect argument
               throw if feature does not have attached slice
  Caller     : general, transform()
  Status     : Stable

=cut
*/
SeqFeature *SeqFeature_transfer(SeqFeature *sf, Slice *slice) {
//  if(!$slice || !ref($slice) || (!$slice->isa('Bio::EnsEMBL::Slice') && !$slice->isa('Bio::EnsEMBL::LRGSlice'))) {
  if (slice == NULL) {
    fprintf(stderr, "Slice argument is required\n");
    exit(1);
  }

  // make a shallow copy of the feature to be transfered
  SeqFeature *feature;
// NIY: How to do the copy, that is the question???
// Temporary hack
  feature = SeqFeatureFactory_newFeatureFromFeature(sf);
  SeqFeature_setStart(feature, SeqFeature_getStart(sf));
  SeqFeature_setEnd(feature, SeqFeature_getEnd(sf));
  SeqFeature_setStrand(feature, SeqFeature_getStrand(sf));
  SeqFeature_setSlice(feature, SeqFeature_getSlice(sf));

  SeqFeature_setDbID(feature, SeqFeature_getDbID(sf));
  //%{$feature} = %{$self};
  //bless $feature, ref($self);
  //weaken $feature->{adaptor};

  Slice *currentSlice = SeqFeature_getSlice(sf);

  if (currentSlice == NULL) {
    fprintf(stderr, "Warning: Feature cannot be transfered without attached slice.\n");
    return NULL;
  }

  CoordSystem *curCs = Slice_getCoordSystem(currentSlice);
  CoordSystem *destCs = Slice_getCoordSystem(slice);

  //if we are not in the same coord system a transformation step is needed first
  if (CoordSystem_compare(destCs, curCs)) {
    SeqFeature *transformedFeature = SeqFeature_transform(feature, CoordSystem_getName(destCs), CoordSystem_getVersion(destCs), slice);
    if (transformedFeature == NULL) {
// NIY: Free feature???
      return NULL;
    }
// NIY: Free feature???
    feature = transformedFeature;
    currentSlice = SeqFeature_getSlice(feature);
  }

  // feature went to entirely different seq_region
  // NIY: Can I use Ids here rather than names????
  if (strcmp(Slice_getSeqRegionName(currentSlice), Slice_getSeqRegionName(slice))) {
// NIY: Anything to free??
    return NULL;
  }

  //if the current feature positions are not relative to the start of the
  //seq region, convert them so they are
  long curSliceStart  = Slice_getStart(currentSlice);
  long curSliceStrand = Slice_getStrand(currentSlice);
  if (curSliceStart != 1 || curSliceStrand != 1) {
    long fStart = SeqFeature_getStart(feature);
    long fEnd   = SeqFeature_getEnd(feature);

    if(curSliceStrand == 1) {
      SeqFeature_setStart(feature, fStart + curSliceStart - 1);
      SeqFeature_setEnd  (feature, fEnd   + curSliceStart - 1);
    } else {
      long curSliceEnd = Slice_getEnd(currentSlice);
      SeqFeature_setStart (feature, curSliceEnd - fEnd   + 1);
      SeqFeature_setEnd   (feature, curSliceEnd - fStart + 1);
      SeqFeature_setStrand(feature, (SeqFeature_getStrand(feature) * -1));
    }
  }

  long fStart = SeqFeature_getStart(feature);
  long fEnd   = SeqFeature_getEnd(feature);

  //convert to destination slice coords
  if (Slice_getStrand(slice) == 1) {
    SeqFeature_setStart (feature, fStart - Slice_getStart(slice) + 1);
    SeqFeature_setEnd   (feature, fEnd   - Slice_getStart(slice) + 1);
  } else {
    SeqFeature_setStart (feature, Slice_getEnd(slice) - fEnd   + 1);
    SeqFeature_setEnd   (feature, Slice_getEnd(slice) - fStart + 1);
    SeqFeature_setStrand(feature, (SeqFeature_getStrand(feature) * -1));
  }

  SeqFeature_setSlice(feature, slice);

  return feature;
}

/*
=head2 project_to_slice

  Arg [1]    : slice to project to


  Example    :
    my $clone_projection = $feature->project_to_slice($slice);

    foreach my $seg (@$clone_projection) {
      my $clone = $seg->to_Slice();
      print "Features current coords ", $seg->from_start, '-',
        $seg->from_end, " project onto clone coords " .
        $clone->seq_region_name, ':', $clone->start, '-', $clone->end,
        $clone->strand, "\n";
    }
  Description: Returns the results of 'projecting' this feature onto another
               slice . This is useful to see where a feature
               would lie in a coordinate system in which it
               crosses a boundary.

               This method returns a reference to a list of
               Bio::EnsEMBL::ProjectionSegment objects.
               ProjectionSegments are blessed arrays and can also be used as
               triplets [from_start,from_end,to_Slice]. The from_start and
               from_end are the coordinates relative to the feature start.
               For example, if a feature is current 100-200bp on a slice
               then the triplets returned might be:
               [1,50,$slice1],
               [51,101,$slice2]

               The to_Slice is a slice spanning the region on the requested
               coordinate system that this feature projected to.

               If the feature projects entirely into a gap then a reference to
               an empty list is returned.

  Returntype : listref of Bio::EnsEMBL::ProjectionSegments
               which can also be used as [$start,$end,$slice] triplets
  Exceptions : slice does not have an adaptor
  Caller     : general
  Status     : At Risk

=cut
*/
Vector *SeqFeature_projectToSlice(SeqFeature *sf, Slice *toSlice) {
  Slice *slice = SeqFeature_getSlice(sf);

  if (slice == NULL) {
    fprintf(stderr, "Warning: Feature cannot be projected without attached slice.\n");
    return Vector_new();
  }


  //get an adaptor from the attached slice because this feature may not yet
  //be stored and may not have its own adaptor
  SliceAdaptor *sa = (SliceAdaptor *)Slice_getAdaptor(slice);

  if (sa == NULL) {
    fprintf(stderr, "Cannot project feature because associated slice does not have an adaptor.\n");
    exit(1);
  }

  int strand = SeqFeature_getStrand(sf) * Slice_getStrand(slice);

  //fetch by feature always gives back forward strand slice:
  Slice *featSlice = SliceAdaptor_fetchByFeature(sa, sf, 0, 0);
  if (strand == -1) {
    Slice *invFeatSlice = Slice_invert(featSlice);

// NIY: I probably want to free this, but nervous about it
//    Slice_free(featSlice);
    featSlice = invFeatSlice;
  }

  Vector *out = Slice_projectToSlice(featSlice, toSlice);

// Again not sure - actually much less sure about this free
//  Slice_free(featSlice);

  return out;
}

/*
=head2 project

  Arg [1]    : string $name
               The name of the coordinate system to project this feature onto
  Arg [2]    : string $version (optional)
               The version of the coordinate system (such as 'NCBI34') to
               project this feature onto
  Example    :
    my $clone_projection = $feature->project('clone');

    foreach my $seg (@$clone_projection) {
      my $clone = $seg->to_Slice();
      print "Features current coords ", $seg->from_start, '-',
        $seg->from_end, " project onto clone coords " .
        $clone->seq_region_name, ':', $clone->start, '-', $clone->end,
        $clone->strand, "\n";
    }
  Description: Returns the results of 'projecting' this feature onto another
               coordinate system.  This is useful to see where a feature
               would lie in a coordinate system in which it
               crosses a boundary.

               This method returns a reference to a list of
               Bio::EnsEMBL::ProjectionSegment objects.
               ProjectionSegments are blessed arrays and can also be used as
               triplets [from_start,from_end,to_Slice]. The from_start and
               from_end are the coordinates relative to the feature start.
               For example, if a feature is current 100-200bp on a slice
               then the triplets returned might be:
               [1,50,$slice1],
               [51,101,$slice2]

               The to_Slice is a slice spanning the region on the requested
               coordinate system that this feature projected to.

               If the feature projects entirely into a gap then a reference to
               an empty list is returned.

  Returntype : listref of Bio::EnsEMBL::ProjectionSegments
               which can also be used as [$start,$end,$slice] triplets
  Exceptions : slice does not have an adaptor
  Caller     : general
  Status     : Stable

=cut
*/
Vector *SeqFeature_project(SeqFeature *sf, char *csName, char *csVersion) {
  Slice *slice = SeqFeature_getSlice(sf);

  if (slice == NULL) {
    fprintf(stderr, "Warning. Feature cannot be projected without attached slice.\n");
    return Vector_new();
  }

  //get an adaptor from the attached slice because this feature may not yet
  //be stored and may not have its own adaptor
  SliceAdaptor *sa = (SliceAdaptor *)Slice_getAdaptor(slice);

  if (sa == NULL) {
    fprintf(stderr, "Cannot project feature because associated slice does not have an adaptor.\n");
    exit(1);
  }

  int strand = SeqFeature_getStrand(sf) * Slice_getStrand(slice);

  //fetch by feature always gives back forward strand slice:
  Slice *featSlice = SliceAdaptor_fetchByFeature(sa, sf, 0, 0);
  if (strand == -1) {
    Slice *invFeatSlice = Slice_invert(featSlice);

// NIY: I probably want to free this, but nervous about it
//    Slice_free(featSlice);
    featSlice = invFeatSlice;
  }

  Vector *out = Slice_project(featSlice, csName, csVersion);

// NIY: Anything to free???

  return out;
}

/*
=head2 seqname

  Arg [1]    : (optional) $seqname
  Example    : $seqname = $feat->seqname();
  Description: Getter/Setter for the name of the sequence that this feature
               is on. Normally you can get away with not setting this value
               and it will default to the name of the slice on which this
               feature is on.  It is useful to set this value on features which
               do not ordinarily sit on features such as ProteinFeatures which
               sit on peptides.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub seqname {
  my $self = shift;
  
  if(@_) {
    $self->{'seqname'} = shift;
  }

  if(!$self->{'seqname'} && $self->slice()) {
    return $self->slice->name();
  }

  return $self->{'seqname'};
}




=head2 display_id

  Arg [1]    : none
  Example    : print $f->display_id();
  Description: This method returns a string that is considered to be
               the 'display' identifier.  It is overridden by subclasses to
               return an appropriate value for objects of that particular 
               class.  If no appropriate display id is available an empty
               string is returned instead.
  Returntype : string
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub display_id {
  my $self = shift;
  return '';
}


=head2 feature_Slice

  Args       : none
  Example    : $slice = $feature->feature_Slice()
  Description: This is a convenience method to return a slice that covers the
               Area of this feature. The feature start will be at 1 on it, and
               it will have the length of this feature.
  Returntype : Bio::EnsEMBL::Slice or undef if this feature has no attached
               Slice.
  Exceptions : warning if Feature does not have attached slice.
  Caller     : web drawing code
  Status     : Stable

=cut

sub feature_Slice {
  my $self = shift;

  my $slice = $self->slice();

  if(!$slice) {
    warning('Cannot obtain Feature_Slice for feature without attached slice');
    return undef;
  }

  if($slice->isa("Bio::EnsEMBL::StrainSlice")){
    return Bio::EnsEMBL::StrainSlice->new
      (-seq_region_name   => $slice->seq_region_name,
       -seq_region_length => $slice->seq_region_length,
       -coord_system      => $slice->coord_system,
       -start             => $self->seq_region_start(),
       -end               => $self->seq_region_end(),
       -strand            => $self->seq_region_strand(),
       -adaptor           => $slice->adaptor(),
       -strain_name       => $slice->strain_name());
  }
  else{
    return Bio::EnsEMBL::Slice->new
      (-seq_region_name   => $slice->seq_region_name,
       -seq_region_length => $slice->seq_region_length,
       -coord_system      => $slice->coord_system,
       -start             => $self->seq_region_start(),
       -end               => $self->seq_region_end(),
       -strand            => $self->seq_region_strand(),
       -adaptor           => $slice->adaptor());
  }
}


=head2 seq_region_name

  Arg [1]    : none
  Example    : print $feature->seq_region_name();
  Description: Gets the name of the seq_region which this feature is on.
               Returns undef if this Feature is not on a slice.
  Returntype : string or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub seq_region_name {
  my $self = shift;
  my $slice = $self->{'slice'};

  return ($slice) ? $slice->seq_region_name() : undef;
}


=head2 seq_region_length

  Arg [1]    : none
  Example    : print $feature->seq_region_length();
  Description: Returns the length of the seq_region which this feature is on 
               Returns undef if this Feature is not on a slice.
  Returntype : int (unsigned) or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut


sub seq_region_length {
  my $self = shift;
  my $slice = $self->{'slice'};

  return ($slice) ? $slice->seq_region_length() : undef;
}

*/

/*
=head2 seq_region_strand

  Arg [1]    : none
  Example    : print $feature->seq_region_strand();
  Description: Returns the strand of the seq_region which this feature is on 
               (i.e. feature_strand * slice_strand)
               Returns undef if this Feature is not on a slice.
  Returntype : 1,0,-1 or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/


int SeqFeature_getSeqRegionStrand(SeqFeature *sf) {
  Slice *slice = SeqFeature_getSlice(sf);

  if (slice) {
    return Slice_getStrand(slice) * SeqFeature_getStrand(sf);
  } else {
// Not sure what to return in this case
    return STRAND_UNDEF;
  }
  //return ($slice) ? $slice->strand() * $self->{'strand'} : undef;
}


/*
=head2 seq_region_start

  Arg [1]    : none
  Example    : print $feature->seq_region_start();
  Description: Convenience method which returns the absolute start of this
               feature on the seq_region, as opposed to the relative (slice) 
               position.

               Returns undef if this feature is not on a slice.
  Returntype : int or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/

int SeqFeature_getSeqRegionStart(SeqFeature *sf) {
  Slice *slice = SeqFeature_getSlice(sf);

  if (slice != NULL) {
    long start = POS_UNDEF;

    if ( Slice_getStrand(slice) == 1 ) {
      if (SeqFeature_getStart(sf) != POS_UNDEF ) {
//        if ($self->start < 0 && $slice->is_circular) {
//          $start = $slice->seq_region_length + $self->start;
//        } else {
        start = Slice_getStart(slice) + SeqFeature_getStart(sf) - 1;
//        }
      }
    } else {
      if (SeqFeature_getEnd(sf) != POS_UNDEF) {
        start = Slice_getEnd(slice) - SeqFeature_getEnd(sf) + 1;
      }
    }

/* Circular shite
    if (start != POS_UNDEF
         && $slice->is_circular()
         && $start > $slice->seq_region_length() ) {
      $start -= $slice->seq_region_length();
    }
*/
    return start;
  }

  return POS_UNDEF;
}


/*
=head2 seq_region_end

  Arg [1]    : none
  Example    : print $feature->seq_region_end();
  Description: Convenience method which returns the absolute end of this
               feature on the seq_region, as opposed to the relative (slice)
               position.

               Returns undef if this feature is not on a slice.
  Returntype : int or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/

int SeqFeature_getSeqRegionEnd(SeqFeature *sf) {
  Slice *slice = SeqFeature_getSlice(sf);

  if (slice != NULL) {
    long end = POS_UNDEF;

    if ( Slice_getStrand(slice) == 1 ) {
      if (SeqFeature_getEnd(sf) != POS_UNDEF) {
        end = Slice_getStart(slice) + SeqFeature_getEnd(sf) - 1;
      }
    } else {
      if (SeqFeature_getStart(sf) != POS_UNDEF) {
        end = Slice_getEnd(slice) - SeqFeature_getStart(sf) + 1;
      }
    }


/* Circular shite
    if (    defined($end)
         && $slice->is_circular()
         && $end > $slice->seq_region_length() ) {
      $end -= $slice->seq_region_length();
    }
*/

    return end;
  }

  return POS_UNDEF;
}


/*
=head2 coord_system_name

  Arg [1]    : none
  Example    : print $feature->coord_system_name()
  Description: Gets the name of the coord_system which this feature is on.
               Returns undef if this Feature is not on a slice.
  Returntype : string or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub coord_system_name {
  my $self = shift;
  my $slice = $self->{'slice'};
  return ($slice) ? $slice->coord_system_name() : undef;
}


=head2 seq

  Args       : none
  Example    : my $dna_sequence = $simple_feature->seq();
  Description: Returns the dna sequence from the attached slice and 
               attached database that overlaps with this feature.
               Returns undef if there is no slice or no database.
               Returns undef if this feature is unstranded (i.e. strand=0).
  Returntype : undef or string
  Exceptions : warning if this feature is not stranded
  Caller     : general
  Status     : Stable

=cut


sub seq {
  my $self = shift;

  if( ! defined $self->{'slice'} ) {
    return undef;
  }

  if(!$self->strand()) {
    warning("Cannot retrieve sequence for unstranded feature.");
    return undef;
  }

  return $self->{'slice'}->subseq($self->start(), $self->end(),
                                  $self->strand());

}




=head2 get_all_alt_locations

  Arg [1]    : Boolean override flag to force the method to return all 
               Features on the reference sequence as well.
               
  Example    : @features = @{$feature->get_all_alt_locations()};
               foreach $f (@features) {
                 print $f->slice->seq_region_name,' ',$f->start, $f->end,"\n";
               }

  Description: Retrieves shallow copies of this feature in its alternate
               locations.  A feature can be considered to have multiple
               locations when it sits on a alternative structural haplotype
               or when it is on a Pseudo Autosomal Region.  Most features will
               just return a reference to an empty list though.
               The features returned by this method will be on a slice which
               covers the entire alternate region.

               Currently this method does not take into account alternate
               locations on the alternate locations (e.g. a reference
               sequence may have multiple alternate haplotypes.  Asking
               for alternate locations of a feature on one of the alternate
               haplotypes will give you back the reference location, but not
               locations on the other alternate haplotypes).

  Returntype : listref of features of the same type of this feature.
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_alt_locations {
  my $self = shift;
  my $return_all = shift || 0;

  my $slice = $self->{'slice'} or return [];
  my $sa = $slice->adaptor() or return [];

  # get slice of entire region
  $slice = $sa->fetch_by_seq_region_id($slice->get_seq_region_id);

  my $axfa = $sa->db->get_AssemblyExceptionFeatureAdaptor();
  my $axfs = $axfa->fetch_all_by_Slice($slice);

  my (@haps, @alt);

  foreach my $axf (@$axfs) {
    if(uc($axf->type()) eq 'HAP') {
      push @haps, $axf;
    } elsif(uc($axf->type()) =~ 'PAR') {
      push @alt, $axf;
    } elsif( $axf->type() eq "PATCH_FIX"){
      push @haps, $axf;
    } elsif( $axf->type() eq "PATCH_FIX REF"){
      push @haps, $axf  if $return_all > 0 ;
    } elsif( $axf->type() eq "HAP REF" ) {
      push @haps, $axf if $return_all > 0 ;
      # do nothing when you are on REF
    } elsif( $axf->type() eq "PATCH_NOVEL"){
      push @haps, $axf;
    }elsif( $axf->type() eq "PATCH_NOVEL REF"){
      push @haps, $axf  if $return_all > 0 ;
    } else {
      warning("Unknown exception feature type ". $axf->type()."- ignoring.");
    }
  }

  # regions surrounding hap are those of interest, not hap itself
  # convert hap alt. exc. features to regions around haps instead
  foreach my $h (@haps) {
    my $haslice = $h->alternate_slice();
    my $hacs    = $haslice->coord_system();

    if($h->start() > 1 && $haslice->start() > 1) {
      my $aslice = $sa->fetch_by_region($hacs->name(),
                                        $haslice->seq_region_name(),
                                        1,
                                        $haslice->start()-1,
                                        $haslice->strand(),
                                        $hacs->version());

      push @alt, Bio::EnsEMBL::AssemblyExceptionFeature->new
        (-start  => 1,
         -end    => $h->start()-1,
         -alternate_slice => $aslice);
    }

    if($h->end() < $slice->seq_region_length() &&
       $haslice->end < $haslice->seq_region_length()) {
      my $aslice = $sa->fetch_by_region($hacs->name(),
                                        $haslice->seq_region_name(),
                                        $haslice->end()+1,
                                        $haslice->seq_region_length(),
                                        $haslice->strand(),
                                        $hacs->version());

      push @alt, Bio::EnsEMBL::AssemblyExceptionFeature->new
        (-start  => $h->end() + 1,
         -end    => $slice->seq_region_length(),
         -alternate_slice => $aslice);
    }
  }


  # check if exception regions contain our feature

  my @features;

  foreach my $axf (@alt) {
    # ignore other region if feature is not entirely on it
    next if($self->seq_region_start() < $axf->start() ||
            $self->seq_region_end()   > $axf->end());

    # quick shallow copy of the feature
    my $f;
    %$f = %$self;
    bless $f, ref($self);

    my $aslice = $axf->alternate_slice();

    # position feature on entire slice of other region
    
    # Cache seq_region_* to prevent contamination when changing feature coordinates.
    my $seq_region_start = $f->seq_region_start();
    my $seq_region_end = $f->seq_region_end();
    
    $f->{'start'}  = $seq_region_start - $axf->start() + $aslice->start();
    $f->{'end'}    = $seq_region_end   - $axf->start() + $aslice->start();
    $f->{'strand'} *= $aslice->strand();

    $f->{'slice'} = $sa->fetch_by_seq_region_id($aslice->get_seq_region_id());

    push @features, $f;
  }

  return \@features;
}


=head2 overlaps

  Arg [1]    : Bio::EnsEMBL::Feature $f
               The other feature you want to check overlap with this feature
               for.
  Description: This method does a range comparison of this features start and
               end and compares it with another features start and end. It will
               return true if these ranges overlap and the features are on the
               same seq_region.
  Returntype : TRUE if features overlap, FALSE if they don't
  Exceptions : warning if features are on different seq_regions
  Caller     : general
  Status     : Stable

=cut

sub overlaps {
  my $self = shift;
  my $f = shift;

  my $sr1_name = $self->seq_region_name;
  my $sr2_name = $f->seq_region_name;

  if ($sr1_name and $sr2_name and ($sr1_name ne $sr2_name)) {
    warning("Bio::EnsEMBL::Feature->overlaps(): features are on different seq regions.");
    return undef;
  }
  
  return ($self->seq_region_end >= $f->seq_region_start and $self->seq_region_start <= $f->seq_region_end);
}


=head2 get_overlapping_Genes

  Description: Get all the genes that overlap this feature.
  Returntype : list ref of Bio::EnsEMBL::Gene
  Caller     : general
  Status     : UnStable

=cut

sub get_overlapping_Genes{
  my $self = shift;

  my $slice = $self->feature_Slice;
  return $slice->get_all_Genes();
}

# query for absolute nearest.
# select x.display_label, g.gene_id, g.seq_region_start, ABS(cast((32921638 - g.seq_region_end) as signed))  as 'dist' from gene g, xref x where g.display_xref_id = x.xref_id and seq_region_id = 27513 order by ABS(cast((32921638 - g.seq_region_end) as signed)) limit 10;

=head2 get_nearest_Gene

  Description: Get the nearest gene to the feature
  Returntype : Bio::EnsEMBL::Gene
  Caller     : general
  Status     : UnStable

=cut

sub get_nearest_Gene {
  my $self = shift;
  my $stranded = shift;
  my $stream = shift;

  my $ga = Bio::EnsEMBL::Registry->get_adaptor($self->adaptor->db->species,"core","Gene");

  return $ga->fetch_nearest_Gene_by_Feature($self, $stranded, $stream);

}

=head2 summary_as_hash

  Example       : $feature_summary = $feature->summary_as_hash();
  Description   : Retrieves a textual summary of this Feature.
                  Should be overidden by subclasses for specific tweaking
  Returns       : hashref of arrays of descriptive strings
  Status        : Intended for internal use
=cut

sub summary_as_hash {
  my $self = shift;
  my %summary;
  $summary{'ID'} = $self->display_id;
  $summary{'start'} = $self->seq_region_start;
  $summary{'end'} = $self->seq_region_end;
  $summary{'strand'} = $self->strand;
  $summary{'seq_region_name'} = $self->seq_region_name;
  return \%summary;
}

=head2 species

  Example     : $feature->species();
  Description : Shortcut to the feature's DBAdaptor and returns its species name 
  Returntype  : String the species name
  Exceptions  : Thrown if there is no attached adaptor
  Caller      : Webcode

=cut

sub species {
  my ($self) = @_;
  throw "Can only call this method if you have attached an adaptor" if ! $self->adaptor();
  return $self->adaptor()->db()->species();
}


##############################################
# Methods included for backwards compatibility
##############################################


=head2 contig

 Deprecated - Included for backwards compatibility only.
 Use slice() instead
=cut
sub contig {
  deprecate('Use slice() instead');
  slice(@_);
}



=head2 sub_SeqFeature

 Deprecated - For genebuild backwards compatibility.
 Avoid using it if possible
=cut
sub sub_SeqFeature{
  my ($self) = @_;
  return @{$self->{'_gsf_sub_array'}} if($self->{'_gsf_sub_array'});
}

=head2 add_sub_SeqFeature

 Deprecated - only for genebuild backward compatibility.
 Avoid using it if possible
=cut
sub add_sub_SeqFeature{
  my ($self,$feat,$expand) = @_;
  my ($p, $f, $l) = caller;
  if( $expand eq 'EXPAND' ) {
    # if this doesn't have start/end set - forget it!
    if( ! $self->start && ! $self->end ) {
      
      $self->start($feat->start());
      $self->end($feat->end());
      $self->strand($feat->strand);
    } else {
      if( $feat->start < $self->start ) {
        $self->start($feat->start);
      }

      if( $feat->end > $self->end ) {
        $self->end($feat->end);
      }
    }
   } else {
     if($self->start > $feat->start || $self->end < $feat->end) {
       throw("$feat is not contained within parent feature, " .
             "and expansion is not valid");
     }
   }

   push(@{$self->{'_gsf_sub_array'}},$feat);
}

=head2 flush_sub_SeqFeature

 Deprecated - Only for genebuild backwards compatibility.
 Avoid using it if possible
=cut
sub flush_sub_SeqFeature {
  my ($self) = @_;
  $self->{'_gsf_sub_array'} = [];
}


sub _deprecated_transform {
  my $self = shift;
  my $arg = shift;

  if(!$arg) {
    warning("Calling transform() with no arguments is deprecated.\n".
          "A coordinate system name argument should be used instead.\n".
          "You probably wanted transform('seqlevel') or transform('contig').");
    return $self->transform('seqlevel');
  }

  if(ref($arg) eq 'Bio::EnsEMBL::Slice') {
    if($arg->{'empty'}) {
      warning("Calling transform with an empty slice is deprecated.\n" .
                "A coordinate system name argument should be used instead.\n".
                "You probably wanted transform('chromosome') or " .
                "transform('toplevel')");
      return $self->transform('toplevel');
    }
    warning("Calling transform with a slice is deprecated.\n" .
              "Use the transfer method instead");
    return $self->transfer($arg);
  }

  warning("Calling transform with a [".ref($arg)."] arg is no longer " .
          "(or never was) supported.  Doing nothing instead.");

  return $self;
}


=head2 id

Deprecated - only included for backwards compatibility.
Use display_id, hseqname, dbID or stable_id instead

=cut

sub id {
  my $self = shift;
  deprecate("id method is not used - use display_id instead");
  return $self->{'stable_id'} if($self->{'stable_id'});
  return $self->{'hseqname'} if($self->{'hseqname'});
  return $self->{'seqname'}  if($self->{'seqname'});
  return $self->{'dbID'};
}

END OF Perl Feature.pm code
*/


Vector *SeqFeature_transformToRawContigImpl(SeqFeature *sf) {
  BaseContig *featContig = SeqFeature_getContig(sf);

  if (featContig) {
    if (featContig->objectType == CLASS_RAWCONTIG) {
      // we are already in rawcontig coords, nothing needs to be done
      Vector_setElementAt(singleEntryVector, 0, sf);
      return singleEntryVector;
    } else if (featContig->objectType == CLASS_SLICE) {
      // transform to raw_contig coords from Slice coords
      return SeqFeature_transformSliceToRawContig(sf);
    } else {
      // Unknown contig type
      fprintf(stderr, "Error: Cannot transform unknown contig type %d\n",featContig->objectType);
      exit(1);
    }
  } else {
    // Can't convert to rawcontig coords without a contig to work with
    fprintf(stderr, "Error: Objects contig is not defined - cannot transform\n");
    exit(1);
  }
}

SeqFeature *SeqFeature_transformToSliceImpl(SeqFeature *sf, Slice *slice) {
  BaseContig *featContig = SeqFeature_getContig(sf);

  if (featContig) {
    if (featContig->objectType == CLASS_RAWCONTIG)  {
      // transform to slice coords from raw contig coords
      return SeqFeature_transformRawContigToSlice(sf, slice);
    } else if (featContig->objectType == CLASS_SLICE) {
      // transform to slice coords from other slice coords
      fprintf(stderr, "Error: Transforms between slices have not been implemented - Nag Steve\n");
      exit(1);
    } else {
      // Unknown contig type
      fprintf(stderr, "Error: Cannot transform unknown contig type %d\n",featContig->objectType);
      exit(1);
    }
  } else {
    // Can't convert to slice coords without a contig to work with
    fprintf(stderr, "Error: Objects contig is not defined - cannot transform\n");
    exit(1);
  }
}

SeqFeature *SeqFeature_transformRawContigToSliceImpl(SeqFeature *sf, Slice *slice) {
  DBAdaptor *dba;
  RawContigAdaptor *rca;
  AssemblyMapperAdaptor *ama;
  AssemblyMapper *assMapper;
  RawContig *rc;
  MapperRangeSet *mapped;
  MapperCoordinate *mc;

  if (!SeqFeature_getContig(sf)) {
    fprintf(stderr, "Error: can't transform coordinates of %p without a contig defined\n", sf);
    exit(1);
  }

  rc = (RawContig *)SeqFeature_getContig(sf);

  if (!RawContig_getAdaptor(rc)) {
    fprintf(stderr, "Error: can't transform coordinates of %p without an adaptor "
                    "attached to the feature's contig\n", sf);
    exit(1);
  }
  dba = RawContig_getAdaptor(rc)->dba;

  ama = DBAdaptor_getAssemblyMapperAdaptor(dba);
  assMapper = AssemblyMapperAdaptor_fetchByType(ama, CoordSystem_getName(Slice_getCoordSystem(slice)));

  rca = DBAdaptor_getRawContigAdaptor(dba);

  mapped = AssemblyMapper_mapCoordinatesToAssembly(
    assMapper,
    //RawContig_getDbID(rc),
    RawContig_getName(rc),
    SeqFeature_getStart(sf),
    SeqFeature_getEnd(sf),
    SeqFeature_getStrand(sf)
  );

  if (!mapped || mapped->nRange == 0) {
    fprintf(stderr, "ERROR: SeqFeature couldnt map dbID = " IDFMTSTR "\n",SeqFeature_getDbID(sf));
    exit(1);
  }

  if (mapped->nRange != 1) {
    fprintf(stderr, "Error: seq feature should only map to one chromosome - "
                    "something bad has happened ...\n");
    exit(1);
  }


  if (MapperRangeSet_getRangeAt(mapped,0)->rangeType == MAPPERRANGE_GAP) {
    fprintf(stderr, "Warning: feature lies on gap\n");
    return NULL;
  }

  mc = (MapperCoordinate *)MapperRangeSet_getRangeAt(mapped,0);


  // the slice is an empty slice, create an enitre chromosome slice and
  // replace the empty slice with it
//  if (Slice_getEmptyFlag(slice)) {
//    SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(dba);
//    ChromosomeAdaptor *ca = DBAdaptor_getChromosomeAdaptor(dba);
//    char *chrName = Chromosome_getName(ChromosomeAdaptor_fetchByDbID(ca,mc->id));
//
//    slice = SliceAdaptor_fetchByChrName(sa, chrName);
//  }

  // mapped coords are on chromosome - need to convert to slice
  if(Slice_getStrand(slice) == 1) {
    SeqFeature_setStart(sf, mc->start - Slice_getChrStart(slice) + 1);
    SeqFeature_setEnd(sf, mc->end   - Slice_getChrStart(slice) + 1);
    SeqFeature_setStrand(sf, mc->strand);
  } else {
    SeqFeature_setStart(sf, Slice_getChrEnd(slice) - mc->end   + 1);
    SeqFeature_setEnd(sf, Slice_getChrEnd(slice) - mc->start + 1);
    SeqFeature_setStrand(sf, mc->strand * -1);
  }

  //set the contig to the slice
  SeqFeature_setContig(sf, slice);

  return sf;
}

Vector *SeqFeature_transformSliceToRawContigImpl(SeqFeature *sf) {
  Slice *slice;
  DBAdaptor *dba;
  RawContigAdaptor *rca;
  AssemblyMapperAdaptor *ama;
  AssemblyMapper *assMapper;
  Vector *out;
  int start;
  int end;
  int strand;
  MapperRangeSet *mapped;
  MapperCoordinate *mc;

  if (!SeqFeature_getContig(sf)) {
    fprintf(stderr, "Error: can't transform coordinates of %p without a contig defined\n", sf);
    exit(1);
  }

  slice = (Slice *)SeqFeature_getContig(sf);

  if (!Slice_getAdaptor(slice)) {
    fprintf(stderr, "Error: can't transform coordinates of %p without an adaptor "
                    "attached to the feature's slice\n", sf);
    exit(1);
  }

  dba = Slice_getAdaptor(slice)->dba;

  ama = DBAdaptor_getAssemblyMapperAdaptor(dba);
  assMapper = AssemblyMapperAdaptor_fetchByType(ama, CoordSystem_getName(Slice_getCoordSystem(slice)));

  rca = DBAdaptor_getRawContigAdaptor(dba);

  // first convert the features coordinates to assembly coordinates
  if (Slice_getStrand(slice) == 1) {
    start  = Slice_getChrStart(slice) + SeqFeature_getStart(sf) - 1;
    end    = Slice_getChrStart(slice) + SeqFeature_getEnd(sf)   - 1;
    strand = SeqFeature_getStrand(sf);
  } else {
    start  = Slice_getChrEnd(slice) - SeqFeature_getEnd(sf)   + 1;
    end    = Slice_getChrEnd(slice) - SeqFeature_getStart(sf) + 1;
    strand = SeqFeature_getStrand(sf) * -1;
  }

  // convert the assembly coordinates to RawContig coordinates
  mapped = AssemblyMapper_mapCoordinatesToRawContig(assMapper,
    Slice_getChrName(slice), //Slice_getChrId(slice),
    start,
    end,
    strand
  );

  if (!mapped || mapped->nRange == 0) {
    fprintf(stderr, "ERROR: SeqFeature couldnt map dbID = " IDFMTSTR "\n",SeqFeature_getDbID(sf));
    exit(1);
  }


  if (mapped->nRange == 1) {
    RawContig *rc;

    if (MapperRangeSet_getRangeAt(mapped,0)->rangeType == MAPPERRANGE_GAP) {
      fprintf(stderr, "Warning: feature lies on gap\n");
      Vector_setElementAt(singleEntryVector,0,sf);
      return singleEntryVector;
    }

    mc = (MapperCoordinate *)MapperRangeSet_getRangeAt(mapped,0);

    rc = RawContigAdaptor_fetchByDbID(rca, mc->id);

    SeqFeature_setStart (sf, mc->start);
    SeqFeature_setEnd   (sf, mc->end);
    SeqFeature_setStrand(sf, mc->strand);
    // NIY SeqFeature_set$self->seqname   ($mapped[0]->id);
    //printf(stderr, "setting contig to be ".$mapped[0]->id."\n";
    SeqFeature_setContig(sf, rc);

    Vector_setElementAt(singleEntryVector,0,sf);
    return singleEntryVector;

  } else {
    Vector *gaps   = Vector_new();
    Vector *coords = Vector_new();
    int i;

    // more than one object returned from mapper
    // possibly more than one RawContig in region

    for (i=0; i<mapped->nRange; i++) {
      MapperRange *mr = MapperRangeSet_getRangeAt(mapped, i);
    
      if (mr->rangeType == MAPPERRANGE_GAP) {
        Vector_addElement(gaps, mr);
      } else if (mr->rangeType == MAPPERRANGE_COORD) {
        Vector_addElement(coords, mr);
      } else {
        fprintf(stderr,"Error: Unknown range type\n");
        exit(1);
      }
    }

    // case where only one RawContig maps
    if (Vector_getNumElement(coords) == 1) {
      mc = (MapperCoordinate *)Vector_getElementAt(coords,0);
 
      SeqFeature_setStart (sf, mc->start);
      SeqFeature_setEnd   (sf, mc->end);
      SeqFeature_setStrand(sf, mc->strand);
      // NIY $self->seqname($coords[0]->id);
      //print STDERR "2 setting contig to be ".$coords[0]->id."\n";
      SeqFeature_setContig(sf, RawContigAdaptor_fetchByDbID(rca, mc->id));

      fprintf(stderr, "Warning: Feature [%p] truncated as lies partially on a gap\n", sf);
      Vector_setElementAt(singleEntryVector,0,sf);
      out = singleEntryVector;

    } else {
      if (SeqFeature_getIsSplittable(sf)) {
        fprintf(stderr, "Warning: Feature spans >1 raw contig - can't split\n");
        Vector_setElementAt(singleEntryVector,0,sf);
// NIY check that this should be a return
        out = singleEntryVector;
      } else {
  
        out = Vector_new();
    
        for (i=0; i<mapped->nRange; i++) {
          SeqFeature *feat;
          MapperCoordinate *mc;
          MapperRange *mr = MapperRangeSet_getRangeAt(mapped,i);
    
          if (mr->rangeType == MAPPERRANGE_GAP) {
            fprintf(stderr, "Warning: piece of seq feature lies on gap\n");
            continue;
          }
    
          mc = (MapperCoordinate *)mr;
  
          feat = SeqFeatureFactory_newFeatureFromFeature(sf);
    
          SeqFeature_setStart(feat, mc->start);
          SeqFeature_setEnd(feat, mc->end);
          SeqFeature_setStrand(feat, mc->strand);
          fprintf(stderr, "3 setting contig to be " IDFMTSTR "\n",mc->id);
          SeqFeature_setContig(feat, RawContigAdaptor_fetchByDbID(rca, mc->id));
          if (SeqFeature_getAdaptor(sf)) SeqFeature_setAdaptor(sf, SeqFeature_getAdaptor(sf));


    // HACK HACK HACK
          if (Class_isDescendent(CLASS_SIMPLEFEATURE, sf->objectType)) {
            SimpleFeature_setDisplayLabel((SimpleFeature *)feat, SimpleFeature_getDisplayLabel((SimpleFeature *)sf));
          }
          SeqFeature_setAnalysis(feat,SeqFeature_getAnalysis(sf));
    
          Vector_addElement(out,feat);
        }
      }
    }
    //NIY freeing mapper and coord etc.
    Vector_free(coords);
    Vector_free(gaps);

    return out;
  }
}

