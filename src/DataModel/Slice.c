#define __SLICE_MAIN__
#include "Slice.h"
#undef __SLICE_MAIN__

#include "DNAAlignFeatureAdaptor.h"
#include "GeneAdaptor.h"
#include "PredictionTranscriptAdaptor.h"
#include "ProteinAlignFeatureAdaptor.h"
#include "RepeatFeatureAdaptor.h"
#include "SimpleFeatureAdaptor.h"
#include "SliceAdaptor.h"
#include "CoordSystemAdaptor.h"
#include "AssemblyMapperAdaptor.h"

#include "Gene.h"

#include "Vector.h"
#include "StrUtil.h"
#include "ProjectionSegment.h"
/*
=head1 DESCRIPTION

A Slice object represents a region of a genome.  It can be used to retrieve
sequence or features from an area of interest.

*/

/*
=head2 new

  Arg [...]  : List of named arguments
               Bio::EnsEMBL::CoordSystem COORD_SYSTEM
               string SEQ_REGION_NAME,
               int    START,
               int    END,
               int    SEQ_REGION_LENGTH, (optional)
               string SEQ (optional)
               int    STRAND, (optional, defaults to 1)
               Bio::EnsEMBL::DBSQL::SliceAdaptor ADAPTOR (optional)
  Example    : $slice = Bio::EnsEMBL::Slice->new(-coord_system => $cs,
                                                 -start => 1,
                                                 -end => 10000,
                                                 -strand => 1,
                                                 -seq_region_name => 'X',
                                                 -seq_region_length => 12e6,
                                                 -adaptor => $slice_adaptor);
  Description: Creates a new slice object.  A slice represents a region
               of sequence in a particular coordinate system.  Slices can be
               used to retrieve sequence and features from an area of
               interest in a genome.

               Coordinates start at 1 and are inclusive.  Negative
               coordinates or coordinates exceeding the length of the
               seq_region are permitted.  Start must be less than or equal.
               to end regardless of the strand.

               Slice objects are immutable. Once instantiated their attributes
               (with the exception of the adaptor) may not be altered.  To
               change the attributes a new slice must be created.
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : throws if start, end, coordsystem or seq_region_name not specified or not of the correct type
  Caller     : general, Bio::EnsEMBL::SliceAdaptor
  Status     : Stable

=cut
*/

Slice *Slice_new(char *regionName, long start, long end, int strand, long length, CoordSystem *coordSystem, SliceAdaptor *sa) {
  Slice *slice;

  if ((slice = (Slice *)calloc(1,sizeof(Slice))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for slice\n");
    return NULL;
  }

  slice->objectType = CLASS_SLICE;
  Object_incRefCount(slice);

  slice->funcs = &sliceFuncs;

/*
  #empty is only for backwards compatibility
  if ($empty) {
    deprecate(   "Creation of empty slices is no longer needed"
               . "and is deprecated" );
     my $self = bless( { 'empty' => 1 }, $class );
    $self->adaptor($adaptor);
    return $self;
  }
*/

  if ( !regionName) {
    fprintf(stderr,"ERROR: Do not have all the parameters for slice\n");
    exit(1);
  }

  if ( start == POS_UNDEF ) { 
    fprintf(stderr, "START argument is POS_UNDEF\n");
    exit(1);
  }
  if ( end   == POS_UNDEF) {
    fprintf(stderr, "END argument is POS_UNDEF\n");
    exit(1);
  }

  if ( start > end + 1 ) {
    fprintf(stderr,"start must be less than or equal to end+1");
    exit(1);
  }

  //if ( !defined($seq_region_length) ) { $seq_region_length = $end }

  if ( length <= 0 ) {
    fprintf(stderr, "SEQ_REGION_LENGTH must be > 0");
    exit(1);
  }

/*
  if ( defined($seq) && CORE::length($seq) != ( $end - $start + 1 ) ) {
    throw('SEQ must be the same length as the defined LENGTH not '
        . CORE::length($seq)
        . ' compared to '
        . ( $end - $start + 1 ) );
  }
*/

  if (coordSystem) {
    if (CoordSystem_getIsTopLevel(coordSystem)) {
      fprintf(stderr, "Cannot create slice on toplevel CoordSystem.");
      exit(1);
    }
  } else {
    fprintf(stderr, "Warning: Slice without coordinate system");
  }

  if (strand == STRAND_UNDEF || strand == 0) strand = 1;

  if (strand != 1 && strand != -1) {
    fprintf(stderr, "STRAND argument must be -1 or 1");
    exit(1);
  }


  Slice_setSeqRegionName(slice,regionName);
  Slice_setSeqRegionStart(slice,start);
  Slice_setSeqRegionEnd(slice,end);
  Slice_setStrand(slice,strand);
  Slice_setSeqRegionLength(slice, length);
  Slice_setCoordSystem(slice, coordSystem);
  Slice_setAdaptor(slice,(BaseAdaptor *)sa);

  //              'seq'               => $seq,
 
  slice->isReference = slice->isTopLevel = slice->hasKaryotype = CHARFLAG_UNSET;

  return slice;
}

ECOSTRING Slice_setSeqRegionName(Slice *sl, char *seqRegionName) {
  EcoString_copyStr(ecoSTable, &(sl->seqRegionName), seqRegionName, 0);

  return sl->seqRegionName;
}


void Slice_free(Slice *slice) {
  Object_decRefCount(slice);

  if (Object_getRefCount(slice) > 0) {
    return;
  } else if (Object_getRefCount(slice) < 0) {
    fprintf(stderr,"Error: Negative reference count for Slice\n"
                   "       Freeing it anyway\n");
  }

  BaseContig_freePtrs(slice);

  if (slice->name) EcoString_freeStr(ecoSTable, slice->name);
  if (slice->seqRegionName) EcoString_freeStr(ecoSTable, slice->seqRegionName);

  free(slice);
}


/*
=head2 coord_system_name

  Arg [1]    : none
  Example    : print $slice->coord_system_name()
  Description: Convenience method.  Gets the name of the coord_system which
               this slice is on.
               Returns undef if this Slice does not have an attached
               CoordSystem.
  Returntype: string or undef
  Exceptions: none
  Caller    : general
  Status     : Stable

=cut
*/

char *Slice_getCoordSystemName(Slice *slice) {
  CoordSystem *cSystem = Slice_getCoordSystem(slice);
  if (cSystem != NULL) {
    return CoordSystem_getName(cSystem);
  } else {
    return NULL;
  }
}


/*
=head2 centrepoint

  Arg [1]    : none
  Example    : $cp = $slice->centrepoint();
  Description: Returns the mid position of this slice relative to the
               start of the sequence region that it was created on.
               Coordinates are inclusive and start at 1.
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/

long Slice_getCentrePoint(Slice *slice) {
  return (Slice_getStart(slice)+Slice_getEnd(slice))/2;
}




/*
=head2 name

  Arg [1]    : none
  Example    : my $results = $cache{$slice->name()};
  Description: Returns the name of this slice. The name is formatted as a colon
               delimited string with the following attributes:
               coord_system:version:seq_region_name:start:end:strand

               Slices with the same name are equivalent and thus the name can
               act as a hash key.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/



ECOSTRING Slice_getName(Slice *slice) {
  char tmpStr[4096];
  CoordSystem *cs = Slice_getCoordSystem(slice);

  if (!slice->name) {
    sprintf(tmpStr,"%s:%s:%s:%ld:%ld:%d",
            ((cs == NULL) ? "" : CoordSystem_getName(cs)),
            ((cs == NULL) ? "" : CoordSystem_getVersion(cs)),
            Slice_getSeqRegionName(slice),
            Slice_getSeqRegionStart(slice),
            Slice_getSeqRegionEnd(slice),
            Slice_getStrand(slice));
    EcoString_copyStr(ecoSTable,&(slice->name),tmpStr,0);
  }
  return slice->name;
}

/*
=head2 is_reference
  Arg        : none
  Example    : my $reference = $slice->is_reference()
  Description: Returns 1 if slice is a reference  slice else 0
  Returntype : int
  Caller     : general
  Status     : At Risk

=cut
*/

int Slice_isReference(Slice *slice) {
  if ( slice->isReference == CHARFLAG_UNSET ) {
    SliceAdaptor *sa = (SliceAdaptor *)Slice_getAdaptor(slice);
 
    if (sa == NULL) {
      fprintf(stderr,"No slice adaptor in Slice_isReference\n");
      exit(1);
    }

    slice->isReference = SliceAdaptor_isReference(sa, Slice_getSeqRegionId(slice) );
  }

  return slice->isReference;
}

/*
=head2 is_toplevel
  Arg        : none
  Example    : my $top = $slice->is_toplevel()
  Description: Returns 1 if slice is a toplevel slice else 0
  Returntype : int
  Caller     : general
  Status     : At Risk

=cut
*/

int Slice_isTopLevel(Slice *slice) {
  if ( slice->isTopLevel == CHARFLAG_UNSET ) {
    SliceAdaptor *sa = (SliceAdaptor *)Slice_getAdaptor(slice);
 
    if (sa == NULL) {
      fprintf(stderr,"No slice adaptor in Slice_isTopLevel\n");
      exit(1);
    }

    slice->isTopLevel = SliceAdaptor_isTopLevel(sa, Slice_getSeqRegionId(slice) );
  }

  return slice->isTopLevel;
}

/*
=head2 has_karyotype
  Arg        : none
  Example    : my $top = $slice->has_karyotype()
  Description: Returns 1 if slice is part of the karyotype else 0
  Returntype : int
  Caller     : general
  Status     : At Risk

=cut
*/
int Slice_hasKaryotype(Slice *slice) {
  if ( slice->hasKaryotype == CHARFLAG_UNSET ) {
    SliceAdaptor *sa = (SliceAdaptor *)Slice_getAdaptor(slice);
 
    if (sa == NULL) {
      fprintf(stderr,"No slice adaptor in Slice_hasKaryotype\n");
      exit(1);
    }

    slice->hasKaryotype = SliceAdaptor_hasKaryotype(sa, Slice_getSeqRegionId(slice) );
  }

  return slice->hasKaryotype;
}

/*
=head2 is_circular
  Arg        : none
  Example    : my $circ = $slice->is_circular()
  Description: Returns 1 if slice is a circular slice else 0
  Returntype : int
  Caller     : general
  Status     : Stable

=cut

sub is_circular {
  my ($self) = @_;
  my $adaptor = $self->adaptor();
  return 0 if ! defined $adaptor;
  if (! exists $self->{'circular'}) {
    my $id = $adaptor->get_seq_region_id($self);
    $self->{circular} = $adaptor->is_circular($id);
  }
  return $self->{circular};
}
*/

/*
=head2 invert

  Arg [1]    : none
  Example    : $inverted_slice = $slice->invert;
  Description: Creates a copy of this slice on the opposite strand and
               returns it.
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/

Slice *Slice_invert(Slice *slice) {
  Slice *inverted = Slice_new(Slice_getSeqRegionName(slice), Slice_getStart(slice), Slice_getEnd(slice), (Slice_getStrand(slice) * -1), 
                              Slice_getSeqRegionLength(slice), Slice_getCoordSystem(slice), (SliceAdaptor *)Slice_getAdaptor(slice));

  
  // reverse complement any attached sequence
  //reverse_comp(\$s{'seq'}) if($s{'seq'});

  return inverted;
}



/*
=head2 seq

  Arg [1]    : none
  Example    : print "SEQUENCE = ", $slice->seq();
  Description: Returns the sequence of the region represented by this
               slice formatted as a string.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/

char *Slice_getSeq(Slice *slice) {
  // special case for in-between (insert) coordinates
  if (Slice_getStart(slice) == Slice_getEnd(slice)+1) {
    char *seq;
    if ((seq = (char *)calloc(1, sizeof(char))) == NULL) {
      fprintf(stderr, "ERROR: Failed allocating sequence string for insert\n");
      exit(1);
    }
    return seq;
  }

  if (Slice_getAdaptor(slice)) {
    SequenceAdaptor *seqa = DBAdaptor_getSequenceAdaptor(Slice_getAdaptor(slice)->dba);
    return SequenceAdaptor_fetchBySliceStartEndStrand(seqa,slice,1,POS_UNDEF,1);
  }

  // no attached sequence, and no db, so just return Ns
  char *seq;
  if ((seq = (char *)calloc(1,sizeof(char))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating sequence string for insert\n");
    exit(1);
  }
  
  SeqUtil_addNs(seq, Slice_getLength(slice));

  return seq;
}


/*
=head2 subseq

  Arg  [1]   : int $startBasePair
               relative to start of slice, which is 1.
  Arg  [2]   : int $endBasePair
               relative to start of slice.
  Arg  [3]   : (optional) int $strand
               The strand of the slice to obtain sequence from. Default
               value is 1.
  Description: returns string of dna sequence
  Returntype : txt
  Exceptions : end should be at least as big as start
               strand must be set
  Caller     : general
  Status     : Stable

=cut
*/

char *Slice_getSubSeq(Slice *slice, int start, int end, int strand) {

  if (end+1 < start) {
    fprintf(stderr,"Error: End coord is less then start coord\n");
    exit(1);
  }

  // special case for in-between (insert) coordinates
  if (end+1 == start) {
    char *seq;
    if ((seq = (char *)calloc(1,sizeof(char))) == NULL) {
      fprintf(stderr, "ERROR: Failed allocating sequence string for insert\n");
      exit(1);
    }
    return seq;
  }

  if (strand != -1 && strand != 1 ) {
    fprintf(stderr,"Error: Invalid strand [%d] in call to Slice subseq.\n",strand);
    exit(1);
  }

  if (Slice_getAdaptor(slice)) {
    SequenceAdaptor *seqa = DBAdaptor_getSequenceAdaptor(Slice_getAdaptor(slice)->dba);
    return SequenceAdaptor_fetchBySliceStartEndStrand(seqa,slice,start,end,strand);
  } else {
    fprintf(stderr, "Not implemented manual seq setting in slice\n");
    exit(1);
/*
    // check for gap at the beginning and pad it with Ns
    if ($start < 1) {
      $subseq = "N" x (1 - $start);
      $start = 1;
    }
    $subseq .= substr ($self->seq(), $start-1, $end - $start + 1);
    // check for gap at the end and pad it with Ns
    if ($end > $self->length()) {
      $subseq .= "N" x ($end - $self->length());
    }
    reverse_comp(\$subseq) if($strand == -1);
*/
  }
}


/*
=head2 sub_Slice_Iterator

  Arg[1]      : int The chunk size to request
  Example     : my $i = $slice->sub_Slice_Iterator(60000); 
                while($i->has_next()) { warn $i->next()->name(); }
  Description : Returns an iterator which batches subslices of this Slice 
                in the requested chunk size
  Returntype  : Bio::EnsEMBL::Utils::Iterator next() will return the next
                 chunk of Slice
  Exceptions  : None

=cut

sub sub_Slice_Iterator {
  my ($self, $chunk_size) = @_;
  throw "Need a chunk size to divide the slice by" if ! $chunk_size;
  my $here = 1;
  my $end = $self->length();
  my $iterator_sub = sub {
    while($here <= $end) {
      my $there = $here + $chunk_size - 1;
      $there = $end if($there > $end); 
      my $slice = $self->sub_Slice($here, $there);
      $here = $there + 1;
      return $slice;
    }
    return;
  };
  return Bio::EnsEMBL::Utils::Iterator->new($iterator_sub);
}
*/

/*
=head2 assembly_exception_type

  Example     : $self->assembly_exception_type(); 
  Description : Returns the type of slice this is. If it is reference then you
                will get 'REF' back. Otherwise you will get the first
                element from C<get_all_AssemblyExceptionFeatures()>. If no
                assembly exception exists you will get an empty string back.
  Returntype  : String
  Exceptions  : None
  Caller      : Public
  Status      : Beta

=cut

sub assembly_exception_type {
  my ($self) = @_;
  my $type = q{};
  if($self->is_reference()) {
    $type = 'REF';
  }
  else {
    my $assembly_exceptions = $self->get_all_AssemblyExceptionFeatures();
    if(@{$assembly_exceptions}) {
      $type = $assembly_exceptions->[0]->type();
    }
  }
  return $type;
}
*/

/*
=head2 is_chromosome

  Example                        : print ($slice->is_chromosome()) ? 'I am a chromosome' : 'Not one'; 
  Description        : Uses a number of rules known to indicate a chromosome region 
                other and takes into account those regions which can be 
                placed on a Chromsome coordinate system but in fact are not
                assembled into one.
  Returntype         : Boolean indicates if the current object is a chromosome
  Exceptions         : None

=cut

sub is_chromosome {
  my ($self) = @_;
  my $coord_system = $self->coord_system->name;
  my $seq_name     = $self->seq_region_name;

  if (($seq_name =~ /random
                    |^Un\d{4}$
                    |^Un\.\d{3}\.\d*$
                    |E\d\d\w*$
                    |_NT_
                    |scaffold_
                    |cutchr
                    |unplaced 
                    |chunk
                    |clone
                    |contig
                    |genescaffold
                    |group
                    |reftig
                    |supercontig
                    |ultracontig        
                    /x) or ( $coord_system !~ /^chromosome$/i )) {
    return 0;
  }
  
  return 1;
}
*/


/*
=head2 get_base_count

  Arg [1]    : none
  Example    : $c_count = $slice->get_base_count->{'c'};
  Description: Retrieves a hashref containing the counts of each bases in the
               sequence spanned by this slice.  The format of the hash is :
               { 'a' => num,
                 'c' => num,
                 't' => num,
                 'g' => num,
                 'n' => num,
                 '%gc' => num }

               All bases which are not in the set [A,a,C,c,T,t,G,g] are
               included in the 'n' count.  The 'n' count could therefore be
               inclusive of ambiguity codes such as 'y'.
               The %gc is the ratio of GC to AT content as in:
               total(GC)/total(ACTG) * 100
               This function is conservative in its memory usage and scales to
               work for entire chromosomes.
  Returntype : hashref
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_base_count {
  my $self = shift;

  my $a = 0;
  my $c = 0;
  my $t = 0;
  my $g = 0;

  my $start = 1;
  my $end;

  my $RANGE = 100_000;
  my $len   = $self->length();

  my $seq;

  while ( $start <= $len ) {
    $end = $start + $RANGE - 1;
    $end = $len if ( $end > $len );

    $seq = $self->subseq( $start, $end );

    $a += $seq =~ tr/Aa//;
    $c += $seq =~ tr/Cc//;
    $t += $seq =~ tr/Tt//;
    $g += $seq =~ tr/Gg//;

    $start = $end + 1;
  }

  my $actg = $a + $c + $t + $g;

  my $gc_content = 0;
  if ( $actg > 0 ) {    # Avoid dividing by 0
    $gc_content = sprintf( "%1.2f", ( ( $g + $c )/$actg )*100 );
  }

  return { 'a'   => $a,
           'c'   => $c,
           't'   => $t,
           'g'   => $g,
           'n'   => $len - $actg,
           '%gc' => $gc_content };
}

*/


/*
=head2 project

  Arg [1]    : string $name
               The name of the coordinate system to project this slice onto
  Arg [2]    : string $version
               The version of the coordinate system (such as 'NCBI34') to
               project this slice onto
  Example    :
    my $clone_projection = $slice->project('clone');

    foreach my $seg (@$clone_projection) {
      my $clone = $segment->to_Slice();
      print $slice->seq_region_name(), ':', $seg->from_start(), '-',
            $seg->from_end(), ' -> ',
            $clone->seq_region_name(), ':', $clone->start(), '-',$clone->end(),
            $clone->strand(), "\n";
    }
  Description: Returns the results of 'projecting' this slice onto another
               coordinate system.  Projecting to a coordinate system that
               the slice is assembled from is analagous to retrieving a tiling
               path.  This method may also be used to 'project up' to a higher
               level coordinate system, however.

               This method returns a listref of triplets [start,end,slice]
               which represents the projection.  The start and end defined the
               region of this slice which is made up of the third value of
               the triplet: a slice in the requested coordinate system.
  Returntype : list reference of Bio::EnsEMBL::ProjectionSegment objects which
               can also be used as [$start,$end,$slice] triplets
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/

Vector *Slice_project(Slice *slice, char *csName, char *csVersion) {

  if (csName == NULL) {
    fprintf(stderr,"Coord_system name argument is required\n");
    exit(1);
  }

  Vector *projection = Vector_new();
  SliceAdaptor *sa = (SliceAdaptor *)Slice_getAdaptor(slice);

  if (sa == NULL) {
    fprintf(stderr, "Warning: Cannot project without attached adaptor.\n");
    return projection;
  }

  if (Slice_getCoordSystem(slice) == NULL) {
    fprintf(stderr, "Warning: Cannot project without attached coord system.\n");
    return projection;
  }


  CoordSystemAdaptor *csa = DBAdaptor_getCoordSystemAdaptor(sa->dba);
  CoordSystem *cs = CoordSystemAdaptor_fetchByName(csa, csName, csVersion);
  CoordSystem *sliceCs = Slice_getCoordSystem(slice);

  if (cs == NULL) {
    fprintf(stderr,"Cannot project to unknown coordinate system [%s %s]", csName, (csVersion != NULL ? csVersion : " undef version"));
    exit(1);
  }

  // no mapping is needed if the requested coord system is the one we are in
  // but we do need to check if some of the slice is outside of defined regions
  if (!CoordSystem_compare(sliceCs, cs)) {
    return Slice_constrainToRegion(slice);
  }

  long currentStart = 1;

  // decompose this slice into its symlinked components.
  // this allows us to handle haplotypes and PARs
  Vector *normalSliceProj = SliceAdaptor_fetchNormalizedSliceProjection(sa, slice, 0);

  int i;
  for (i=0; i<Vector_getNumElement(normalSliceProj); i++) {
    ProjectionSegment *segment = Vector_getElementAt(normalSliceProj, i);

    Slice *normalSlice = ProjectionSegment_getToSlice(segment);

    sliceCs = Slice_getCoordSystem(normalSlice);

    AssemblyMapperAdaptor *asma = DBAdaptor_getAssemblyMapperAdaptor(sa->dba);
    AssemblyMapper *asmMapper = AssemblyMapperAdaptor_fetchByCoordSystems(asma, sliceCs, cs);

    // perform the mapping between this slice and the requested system
    MapperRangeSet *coords;

    if (asmMapper != NULL ) {
      coords = AssemblyMapper_map(asmMapper, Slice_getSeqRegionName(normalSlice), 
                                             Slice_getStart(normalSlice), 
                                             Slice_getEnd(normalSlice), 
                                             Slice_getStrand(normalSlice), 
                                             sliceCs, 0, NULL);
    } else {
      coords = MapperRangeSet_new();
      MapperRangeSet_addRange(coords, (MapperRange *)MapperGap_new(Slice_getStart(normalSlice), Slice_getEnd(normalSlice),0));
    }


    // my $last_rank = 0;
    //construct a projection from the mapping results and return it
    int j;
    for (j=0; j<MapperRangeSet_getNumRange(coords); j++) {
      MapperRange *coord = MapperRangeSet_getRangeAt(coords, j);

      long coordStart  = coord->start;
      long coordEnd    = coord->end;
      long length      = coordEnd - coordStart + 1;

      if (coordStart > coordEnd ) {
        length = Slice_getSeqRegionLength(normalSlice) -
                 coordStart +
                 coordEnd + 1;
      }

//      if( $last_rank != $coord->rank){
//        $current_start = 1;
//        print "LAST rank has changed to ".$coord->rank."from $last_rank \n";
//     }
//      $last_rank = $coord->rank;

      //skip gaps
      if (coord->rangeType == MAPPERRANGE_COORD) {
        MapperCoordinate *mc = (MapperCoordinate *)coord;

        CoordSystem *coordCs = mc->coordSystem;

        // If the normalised projection just ended up mapping to the
        // same coordinate system we were already in then we should just
        // return the original region.  This can happen for example, if we
        // were on a PAR region on Y which refered to X and a projection to
        // 'toplevel' was requested.
        if (!CoordSystem_compare(coordCs,sliceCs)) {
          // trim off regions which are not defined
          return Slice_constrainToRegion(slice);
        }
        //create slices for the mapped-to coord system
        Slice *mcSlice = SliceAdaptor_fetchBySeqRegionId(sa, mc->id, coordStart, coordEnd, mc->strand);

        long currentEnd = currentStart + length - 1;

/*
        if ($current_end > $slice->seq_region_length() && $slice->is_circular ) {
            $current_end -= $slice->seq_region_length();
        }
*/
        Vector_addElement(projection, ProjectionSegment_new(currentStart, currentEnd, mcSlice));
      }

      currentStart += length;
    }
  }

  return projection;
}


Vector *Slice_constrainToRegion(Slice *slice) {

  long entireLen = Slice_getSeqRegionLength(slice);

  // if the slice has negative coordinates or coordinates exceeding the
  // length of the sequence region we want to shrink the slice to
  // the defined region
  Vector *projection = Vector_new();

  if (Slice_getStart(slice) > entireLen || Slice_getEnd(slice) < 1) {
    //none of this slice is in a defined region
    return projection;
  }

  long rightContract = 0;
  long leftContract  = 0;
  if (Slice_getEnd(slice) > entireLen) {
    rightContract = entireLen - Slice_getEnd(slice);
  }
  if (Slice_getStart(slice) < 1) {
    leftContract = Slice_getStart(slice) - 1;
  }

  Slice *newSlice;
  long junk;
  if (leftContract || rightContract) {
    //if slice in negative strand, need to swap contracts
    if (Slice_getStrand(slice) == 1) {
      newSlice = Slice_expand(slice, leftContract, rightContract, 0, &junk, &junk);
    } else if (Slice_getStrand(slice) == -1) {
      newSlice = Slice_expand(slice, rightContract, leftContract, 0, &junk, &junk);
    }
  } else {
    newSlice = slice;
  }

  Vector_addElement(projection, ProjectionSegment_new(1-leftContract, Slice_getLength(slice)+rightContract, newSlice));

  return projection;
}

/*

=head2 expand

  Arg [1]    : (optional) int $five_prime_expand
               The number of basepairs to shift this slices five_prime
               coordinate by.  Positive values make the slice larger,
               negative make the slice smaller.
               coordinate left.
               Default = 0.
  Arg [2]    : (optional) int $three_prime_expand
               The number of basepairs to shift this slices three_prime
               coordinate by. Positive values make the slice larger,
               negative make the slice smaller.
               Default = 0.
  Arg [3]    : (optional) bool $force_expand
               if set to 1, then the slice will be contracted even in the case
               when shifts $five_prime_expand and $three_prime_expand overlap.
               In that case $five_prime_expand and $three_prime_expand will be set
               to a maximum possible number and that will result in the slice
               which would have only 2pbs.
               Default = 0.
  Arg [4]    : (optional) int* $fpref
               The reference to a number of basepairs to shift this slices five_prime
               coordinate by. Normally it would be set to $five_prime_expand.
               But in case when $five_prime_expand shift can not be applied and
               $force_expand is set to 1, then $$fpref will contain the maximum possible
               shift
  Arg [5]    : (optional) int* $tpref
               The reference to a number of basepairs to shift this slices three_prime
               coordinate by. Normally it would be set to $three_prime_expand.
               But in case when $five_prime_expand shift can not be applied and
               $force_expand is set to 1, then $$tpref will contain the maximum possible
               shift
  Example    : my $expanded_slice      = $slice->expand( 1000, 1000);
               my $contracted_slice    = $slice->expand(-1000,-1000);
               my $shifted_right_slice = $slice->expand(-1000, 1000);
               my $shifted_left_slice  = $slice->expand( 1000,-1000);
               my $forced_contracted_slice    = $slice->expand(-1000,-1000, 1, \$five_prime_shift, \$three_prime_shift);

  Description: Returns a slice which is a resized copy of this slice.  The
               start and end are moved outwards from the center of the slice
               if positive values are provided and moved inwards if negative
               values are provided. This slice remains unchanged.  A slice
               may not be contracted below 1bp but may grow to be arbitrarily
               large.
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : warning if an attempt is made to contract the slice below 1bp
  Caller     : general
  Status     : Stable

=cut
*/

Slice *Slice_expand(Slice *slice, long fivePrimeShift, long threePrimeShift, int forceExpand, long *fpRef, long *tpRef) {
/*
  my $self              = shift;
  my $five_prime_shift  = shift || 0;
  my $three_prime_shift = shift || 0;
  my $force_expand      = shift || 0;
  my $fpref             = shift;
  my $tpref             = shift;
*/

/* NIY: Not handling manually added seq yet
  if ( $self->{'seq'} ) {
    warning(
       "Cannot expand a slice which has a manually attached sequence ");
    return undef;
  }
*/

  long sShift = fivePrimeShift;
  long eShift = threePrimeShift;

  if ( Slice_getStrand(slice) != 1 ) {
    eShift = fivePrimeShift;
    sShift = threePrimeShift;
  }

  long newStart = Slice_getStart(slice) - sShift;
  long newEnd   = Slice_getEnd(slice) + eShift;

  /*
  # Wrap around on circular slices
  if (( $new_start <= 0 || $new_start > $self->seq_region_length() || $new_end <= 0 
        || $new_end > $self->seq_region_length() ) && ( $self->is_circular() ) ) {
      
      if ( $new_start <= 0 ) {
        $new_start = $self->seq_region_length() + $new_start;
      }
      if ( $new_start > $self->seq_region_length() ) {
        $new_start -= $self->seq_region_length();
      }
  
      if ( $new_end <= 0 ) {
        $new_end = $self->seq_region_length() + $new_end;
      }
      if ( $new_end > $self->seq_region_length() ) {
        $new_end -= $self->seq_region_length();
      }      
      
  }
  */

  if ( newStart > newEnd) { //  && (not $self->is_circular() ) ) {
    if (forceExpand) {
      // Apply max possible shift, if force_expand is set
      if ( sShift < 0 ) {
        // if we are contracting the slice from the start - move the
        // start just before the end
        newStart = newEnd - 1;
        sShift    = Slice_getStart(slice) - newStart;
      }

      if ( newStart > newEnd ) {
        // if the slice still has a negative length - try to move the
        // end
        if ( eShift < 0 ) {
          newEnd = newStart + 1;
          eShift = newEnd - Slice_getEnd(slice);
        }
      }

      // return the values by which the primes were actually shifted
      *tpRef = (Slice_getStrand(slice) == 1 ? eShift : sShift);
      *fpRef = (Slice_getStrand(slice) == 1 ? sShift : eShift);
    }
    if ( newStart > newEnd ) {
      fprintf(stderr, "Slice start cannot be greater than slice end\n");
    }
  }

  Slice *newSlice = Slice_new(Slice_getSeqRegionName(slice), newStart, newEnd, Slice_getStrand(slice),
                              Slice_getSeqRegionLength(slice), Slice_getCoordSystem(slice), (SliceAdaptor *)Slice_getAdaptor(slice));

  return newSlice;
}

/*
=head2 constrain_to_seq_region
  Example    : $new_slice = $slice->expand(1000,10000);
               $new_slice = $new_slice->constrain_to_seq_region();
  Description: Used to prevent overly zealous expand calls going off the end of
               the sequence region. It contracts the start and end where needed
               and produces a slice copy with the tweaked coordinates.
  Returntype : Bio::EnsEMBL::Slice
=cut
*/

Slice *Slice_constrainToSeqRegion(Slice *slice) {
    // circular calculations should already be taken care of
    //if ($self->is_circular) {return $self} 

  long newStart = Slice_getStart(slice);
  long newEnd   = Slice_getEnd(slice);
    
/* This seems odd -seqRegionSlice returns a 1 to seqRegionLength slice so why can't I just use those values to check against rather than making this unnecessary slice
   Maybe in one of the other obscure slice types, but I'm not going to care about it for now
    my $seq_region = $self->seq_region_Slice;
    
    if ($new_start < $seq_region->start) {$new_start = $seq_region->start}
    if ($new_end > $seq_region->end) {$new_end = $seq_region->end}
*/

  if (newStart < 1) newStart = 1;
  if (newEnd > Slice_getSeqRegionLength(slice)) newEnd = Slice_getSeqRegionLength(slice);
    
    
  Slice *newSlice = Slice_new(Slice_getSeqRegionName(slice), newStart, newEnd, Slice_getStrand(slice),
                              Slice_getSeqRegionLength(slice), Slice_getCoordSystem(slice), (SliceAdaptor *)Slice_getAdaptor(slice));

  return newSlice;
}


/*
=head2 sub_Slice

  Arg   1    : int $start
  Arg   2    : int $end
  Arge [3]   : int $strand
  Example    : none
  Description: Makes another Slice that covers only part of this Slice
               If a Slice is requested which lies outside of the boundaries
               of this function will return undef.  This means that
               behaviour will be consistant whether or not the slice is
               attached to the database (i.e. if there is attached sequence
               to the slice).  Alternatively the expand() method or the
               SliceAdaptor::fetch_by_region method can be used instead.
  Returntype : Bio::EnsEMBL::Slice or undef if arguments are wrong
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/

// Should use 1 as default strand if called without the arg
Slice *Slice_getSubSlice(Slice *slice, long start, long end, int strand) {

  if (start < 1 || start > Slice_getEnd(slice)) {
    //# throw( "start argument not valid" );
    return NULL;
  }

  if( end < start || end > Slice_getEnd(slice)) {
    //# throw( "end argument not valid" )
    return NULL;
  }

  long newStart;
  long newEnd;
  int  newStrand;
  // char *newSeq; // not doing seq yet

  if (Slice_getStrand(slice) == 1 ) {
    newStart = Slice_getStart(slice) + start - 1;
    newEnd = Slice_getStart(slice) + end - 1;
    newStrand = strand;
  } else {
    newStart = Slice_getEnd(slice) - end + 1;;
    newEnd = Slice_getEnd(slice) - start + 1;
    newStrand = -strand;
  }

/*
  if( defined $self->{'seq'} ) {
    $new_seq = $self->subseq( $start, $end, $strand );
  }
*/

  Slice *newSlice = Slice_new(Slice_getSeqRegionName(slice), newStart, newEnd, newStrand,
                              Slice_getSeqRegionLength(slice), Slice_getCoordSystem(slice), (SliceAdaptor *)Slice_getAdaptor(slice));
/*
  if( $new_seq ) {
    $new_slice->{'seq'} = $new_seq;
  }
*/

  return newSlice;
}



/*
=head2 seq_region_Slice

  Arg [1]    : none
  Example    : $slice = $slice->seq_region_Slice();
  Description: Returns a slice which spans the whole seq_region which this slice
               is on.  For example if this is a slice which spans a small region
               of chromosome X, this method will return a slice which covers the
               entire chromosome X. The returned slice will always have strand
               of 1 and start of 1.  This method cannot be used if the sequence
               of the slice has been set manually.
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : warning if called when sequence of Slice has been set manually.
  Caller     : general
  Status     : Stable

=cut
*/

Slice *Slice_getSeqRegionSlice(Slice *slice) {

/*
  if($self->{'seq'}){
    warning("Cannot get a seq_region_Slice of a slice which has manually ".
            "attached sequence ");
    return undef;
  }
*/

  Slice *newSlice = Slice_new(Slice_getSeqRegionName(slice), 1, Slice_getSeqRegionLength(slice), 1,
                              Slice_getSeqRegionLength(slice), Slice_getCoordSystem(slice), (SliceAdaptor *)Slice_getAdaptor(slice));

  return newSlice;
}


/*
=head2 get_seq_region_id

  Arg [1]    : none
  Example    : my $seq_region_id = $slice->get_seq_region_id();
  Description: Gets the internal identifier of the seq_region that this slice
               is on. Note that this function will not work correctly if this
               slice does not have an attached adaptor. Also note that it may
               be better to go through the SliceAdaptor::get_seq_region_id
               method if you are working with multiple databases since is
               possible to work with slices from databases with different
               internal seq_region identifiers.
  Returntype : int or undef if slices does not have attached adaptor
  Exceptions : warning if slice is not associated with a SliceAdaptor
  Caller     : assembly loading scripts, general
  Status     : Stable

=cut
*/

IDType Slice_getSeqRegionId(Slice *slice) {

  if (Slice_getAdaptor(slice)) {
    return SliceAdaptor_getSeqRegionId((SliceAdaptor *)Slice_getAdaptor(slice), slice);
  } else {
    fprintf(stderr, "Warning: Cannot retrieve seq_region_id without attached adaptor.\n");
    return 0;
  }
}


/*
=head2 get_all_Attributes

  Arg [1]    : optional string $attrib_code
               The code of the attribute type to retrieve values for.
  Example    : ($htg_phase) = @{$slice->get_all_Attributes('htg_phase')};
               @slice_attributes    = @{$slice->get_all_Attributes()};
  Description: Gets a list of Attributes of this slice''s seq_region.
               Optionally just get Attrubutes for given code.
  Returntype : listref Bio::EnsEMBL::Attribute
  Exceptions : warning if slice does not have attached adaptor
  Caller     : general
  Status     : Stable

=cut
*/

/* NIY: Need to implement AttributeAdaptor before enabling this
Vector *Slice_getAllAttributes(Slice *slice, char *attribCode) {
  SliceAdaptor *sa = Slice_getAdaptor(slice);
  if (sa == NULL) {
    fprintf(stderr,"Warning: Cannot get attributes without an adaptor.\n");
    return Vector_new();
  }


  AttributeAdaptor *aa = DBAdaptor_getAttributeAdaptor(sa->dba);

  AttributeAdaptor
  Vector *results = AttributeAdaptor_fetchAllBySlice( slice );

  
  if (attribCode!=NULL) {
    char *ucAttribCode;
    StrUtil_copyString(&ucAttribCode, attribCode, 0);
    StrUtil_strupr(ucAttribCode);

    int i;
    for  (i=0;i<Vector_getNumElement(results);i++) {
      Attribute *attrib = Vector_getElementAt(results, i);
      char *ucCode;
      StrUtil_copyString(&ucCode, Attribute_getCode(attrib), 0);
      StrUtil_strupr(ucCode);
      
      if (strcmp(ucCode, ucAttribCode)) {
        Vector_removeElementAt(results, i);
        i--;
      }
      free(ucCode);
    }
    free(ucAttribCode);
  }

  return results;
}
*/


/*
=head2 get_all_PredictionTranscripts

  Arg [1]    : (optional) string $logic_name
               The name of the analysis used to generate the prediction
               transcripts obtained.
  Arg [2]    : (optional) boolean $load_exons
               If set to true will force loading of all PredictionExons
               immediately rather than loading them on demand later.  This
               is faster if there are a large number of PredictionTranscripts
               and the exons will be used.
  Example    : @transcripts = @{$slice->get_all_PredictionTranscripts};
  Description: Retrieves the list of prediction transcripts which overlap
               this slice with logic_name $logic_name.  If logic_name is
               not defined then all prediction transcripts are retrieved.
  Returntype : listref of Bio::EnsEMBL::PredictionTranscript
  Exceptions : warning if slice does not have attached adaptor
  Caller     : none
  Status     : Stable

=cut
*/

DBAdaptor *Slice_getSelectedDBAdaptor(Slice *slice, char *dbType) {
  SliceAdaptor *sa = (SliceAdaptor *)Slice_getAdaptor(slice);

  if (sa == NULL) {
    fprintf(stderr, "Can't get a db adaptor without an adaptor on slice\n");
    exit(1);
  }

  if (dbType) {
    fprintf(stderr, "Selecting dbType not implemented\n");
    return NULL;
  }
  return sa->dba;
}
/*
Vector *Slice_getAllPredictionTranscripts(Slice *slice, char *logicName, int loadExons, char *dbType) {
  if (!Slice_getAdaptor(slice)) {
    fprintf(stderr, "Warning: Cannot get PredictionTranscripts without attached adaptor");
    return Vector_new();
  }

  DBAdaptor *db = Slice_getSelectedDBAdaptor(slice, dbType);

  if (!db) {
    fprintf(stderr, "Warning: Don't have db %s returning empty list\n", dbType);
    return Vector_new();
  }

  PredictionTranscriptAdaptor *pta = DBAdaptor_getPredictionTranscriptAdaptor(db);

  return PredictionTranscriptAdaptor_fetchAllBySlice(pta, slice, logicName, loadExons);
}
*/



/*
=head2 get_all_DnaAlignFeatures

  Arg [1]    : (optional) string $logic_name
               The name of the analysis performed on the dna align features
               to obtain.
  Arg [2]    : (optional) float $score
               The mimimum score of the features to retrieve
  Arg [3]    : (optional) string $dbtype
               The name of an attached database to retrieve the features from
               instead, e.g. 'otherfeatures'.
  Arg [4]    : (optional) float hcoverage
               The minimum hcoverage od the featurs to retrieve
  Example    : @dna_dna_align_feats = @{$slice->get_all_DnaAlignFeatures};
  Description: Retrieves the DnaDnaAlignFeatures which overlap this slice with
               logic name $logic_name and with score above $score.  If
               $logic_name is not defined features of all logic names are
               retrieved.  If $score is not defined features of all scores are
               retrieved.
  Returntype : listref of Bio::EnsEMBL::DnaDnaAlignFeatures
  Exceptions : warning if slice does not have attached adaptor
  Caller     : general
  Status     : Stable

=cut
*/

Vector *Slice_getAllDNAAlignFeatures(Slice *slice, char *logicName, double *score, char *dbType, double *hCoverage) {
  if (!Slice_getAdaptor(slice)) {
    fprintf(stderr, "Warning: Cannot get DnaAlignFeatures without attached adaptor");
    return Vector_new();
  }

  DBAdaptor *db = Slice_getSelectedDBAdaptor(slice, dbType);

  if (!db) {
    fprintf(stderr, "Warning: Don't have db %s returning empty list\n", dbType);
    return Vector_new();
  }

  DNAAlignFeatureAdaptor *dafa = DBAdaptor_getDNAAlignFeatureAdaptor(db);

  if (score && hCoverage){
    fprintf(stderr, "Warning: Can not specify score and hcoverage. Using score only\n");
  }
  if (score){
    return DNAAlignFeatureAdaptor_fetchAllBySliceAndScore(dafa, slice, score, logicName);
  }
    return DNAAlignFeatureAdaptor_fetchAllBySliceAndScore(dafa, slice, score, logicName);
//  return DNAAlignFeatureAdaptor_fetchAllBySliceAndHCoverage(dafa, slice, hCoverage, logicName);
}

/*

=head2 get_all_ProteinAlignFeatures

  Arg [1]    : (optional) string $logic_name
               The name of the analysis performed on the protein align features
               to obtain.
  Arg [2]    : (optional) float $score
               The mimimum score of the features to retrieve
  Arg [3]    : (optional) string $dbtype
               The name of an attached database to retrieve features from
               instead.
  Arg [4]    : (optional) float hcoverage
               The minimum hcoverage od the featurs to retrieve
  Example    : @dna_pep_align_feats = @{$slice->get_all_ProteinAlignFeatures};
  Description: Retrieves the DnaPepAlignFeatures which overlap this slice with
               logic name $logic_name and with score above $score.  If
               $logic_name is not defined features of all logic names are
               retrieved.  If $score is not defined features of all scores are
               retrieved.
  Returntype : listref of Bio::EnsEMBL::DnaPepAlignFeatures
  Exceptions : warning if slice does not have attached adaptor
  Caller     : general
  Status     : Stable

=cut
*/
Vector *Slice_getAllProteinAlignFeatures(Slice *slice, char *logicName, double *score, char *dbType, double *hCoverage) {
  if (!Slice_getAdaptor(slice)) {
    fprintf(stderr, "Warning: Cannot get ProteinAlignFeatures without attached adaptor");
    return Vector_new();
  }

  DBAdaptor *db = Slice_getSelectedDBAdaptor(slice, dbType);

  if (!db) {
    fprintf(stderr, "Warning: Don't have db %s returning empty list\n", dbType);
    return Vector_new();
  }

  ProteinAlignFeatureAdaptor *pafa = DBAdaptor_getProteinAlignFeatureAdaptor(db);

  if (score && hCoverage){
    fprintf(stderr, "Warning: Can not specify score and hcoverage. Using score only\n");
  }
  if (score){
    return ProteinAlignFeatureAdaptor_fetchAllBySliceAndScore(pafa, slice, score, logicName);
  }
    return ProteinAlignFeatureAdaptor_fetchAllBySliceAndScore(pafa, slice, score, logicName);
//  return ProteinAlignFeatureAdaptor_fetchAllBySliceAndHCoverage(pafa, slice, hCoverage, logicName);
}


/*
=head2 get_all_SimilarityFeatures

  Arg [1]    : (optional) string $logic_name
               the name of the analysis performed on the features to retrieve
  Arg [2]    : (optional) float $score
               the lower bound of the score of the features to be retrieved
  Example    : @feats = @{$slice->get_all_SimilarityFeatures};
  Description: Retrieves all dna_align_features and protein_align_features
               with analysis named $logic_name and with score above $score.
               It is probably faster to use get_all_ProteinAlignFeatures or
               get_all_DnaAlignFeatures if a sepcific feature type is desired.
               If $logic_name is not defined features of all logic names are
               retrieved.  If $score is not defined features of all scores are
               retrieved.
  Returntype : listref of Bio::EnsEMBL::BaseAlignFeatures
  Exceptions : warning if slice does not have attached adaptor
  Caller     : general
  Status     : Stable

=cut
*/

/*
Vector *Slice_getAllSimilarityFeatures(Slice *slice, char *logicName, double *score) {

  Vector *out = Slice_getAllProteinAlignFeatures(slice, logicName, score, NULL, NULL);
  Vector *tmp = Slice_getAllDNAAlignFeatures(slice, logicName, score, NULL, NULL);
  Vector_append(out, tmp);
  Vector_free(tmp);

  return out;
}
*/

/*
=head2 get_all_SimpleFeatures

  Arg [1]    : (optional) string $logic_name
               The name of the analysis performed on the simple features
               to obtain.
  Arg [2]    : (optional) float $score
               The mimimum score of the features to retrieve
  Example    : @simple_feats = @{$slice->get_all_SimpleFeatures};
  Description: Retrieves the SimpleFeatures which overlap this slice with
               logic name $logic_name and with score above $score.  If
               $logic_name is not defined features of all logic names are
               retrieved.  If $score is not defined features of all scores are
               retrieved.
  Returntype : listref of Bio::EnsEMBL::SimpleFeatures
  Exceptions : warning if slice does not have attached adaptor
  Caller     : general
  Status     : Stable

=cut
*/

Vector *Slice_getAllSimpleFeatures(Slice *slice, char *logicName, double *score, char *dbType) {
  if (!Slice_getAdaptor(slice)) {
    fprintf(stderr, "Warning: Cannot get SimpleFeatures without attached adaptor");
    return Vector_new();
  }

  DBAdaptor *db = Slice_getSelectedDBAdaptor(slice, dbType);

  if (!db) {
    fprintf(stderr, "Warning: Don't have db %s returning empty list\n", dbType);
    return Vector_new();
  }

  SimpleFeatureAdaptor *sfa = DBAdaptor_getSimpleFeatureAdaptor(db);

  return SimpleFeatureAdaptor_fetchAllBySliceAndScore(sfa, slice, score, logicName);
}


/*
=head2 get_all_RepeatFeatures

  Arg [1]    : (optional) string $logic_name
               The name of the analysis performed on the repeat features
               to obtain.
  Arg [2]    : (optional) string/array $repeat_type
               Limits features returned to those of the specified 
               repeat_type. Can specify a single value or an array reference
               to limit by more than one
  Arg [3]    : (optional) string $db
               Key for database e.g. core/vega/cdna/....
  Example    : @repeat_feats = @{$slice->get_all_RepeatFeatures(undef,'Type II Transposons')};
  Description: Retrieves the RepeatFeatures which overlap  with
               logic name $logic_name and with score above $score.  If
               $logic_name is not defined features of all logic names are
               retrieved.
  Returntype : listref of Bio::EnsEMBL::RepeatFeatures
  Exceptions : warning if slice does not have attached adaptor
  Caller     : general
  Status     : Stable

=cut
*/
Vector *Slice_getAllRepeatFeatures(Slice *slice, char *logicName, char *repeatType, char *dbType) {
  if (!Slice_getAdaptor(slice)) {
    fprintf(stderr, "Warning: Cannot get RepeatFeatures without attached adaptor");
    return Vector_new();
  }

  DBAdaptor *db = Slice_getSelectedDBAdaptor(slice, dbType);

  if (!db) {
    fprintf(stderr, "Warning: Don't have db %s returning empty list\n", dbType);
    return Vector_new();
  }

  RepeatFeatureAdaptor *rpfa = DBAdaptor_getRepeatFeatureAdaptor(db);

  return RepeatFeatureAdaptor_fetchAllBySlice(rpfa, slice, logicName, repeatType);
}

/*
=head2 get_all_LD_values

    Arg [1]     : (optional) Bio::EnsEMBL::Variation::Population $population
    Description : returns all LD values on this slice. This function will only work correctly if the variation
                  database has been attached to the core database. If the argument is passed, will return the LD information
                  in that population
    ReturnType  : Bio::EnsEMBL::Variation::LDFeatureContainer
    Exceptions  : none
    Caller      : contigview, snpview
     Status     : Stable

=cut

sub get_all_LD_values{
    my $self = shift;
    my $population = shift;


    if(!$self->adaptor()) {
        warning('Cannot get LDFeatureContainer without attached adaptor');
        return [];
    }

    my $variation_db = $self->adaptor->db->get_db_adaptor('variation');

    unless($variation_db) {
        warning("Variation database must be attached to core database to " .
                "retrieve variation information" );
        return [];
    }

    my $ld_adaptor = $variation_db->get_LDFeatureContainerAdaptor;

    if( $ld_adaptor ) {
        return $ld_adaptor->fetch_by_Slice($self,$population);
    } else {
        return [];

    }

#     my $ld_adaptor = Bio::EnsEMBL::DBSQL::MergedAdaptor->new(-species => $self->adaptor()->db()->species, -type => "LDFeatureContainer");

#   if( $ld_adaptor ) {
#       my $ld_values = $ld_adaptor->fetch_by_Slice($self,$population);
#       if (@{$ld_values} > 1){
#           warning("More than 1 variation database attached. Trying to merge LD results");
#           my $ld_value_merged = shift @{$ld_values};
#           #with more than 1 variation database attached, will try to merge in one single LDContainer object.
#           foreach my $ld (@{$ld_values}){
#               #copy the ld values to the result hash
#               foreach my $key (keys %{$ld->{'ldContainer'}}){
#                   $ld_value_merged->{'ldContainer'}->{$key} = $ld->{'ldContainer'}->{$key};
#               }
#               #and copy the variationFeatures as well
#               foreach my $key (keys %{$ld->{'variationFeatures'}}){
#                   $ld_value_merged->{'variationFeatures'}->{$key} = $ld->{'variationFeatures'}->{$key};
#               }

#           }
#           return $ld_value_merged;
#       }
#       else{
#           return shift @{$ld_values};
#       }
# } else {
#     warning("Variation database must be attached to core database to " .
#                 "retrieve variation information" );
#     return [];
# }
}

sub _get_VariationFeatureAdaptor {
    
  my $self = shift;
    
  if(!$self->adaptor()) {
    warning('Cannot get variation features without attached adaptor');
    return undef;
  }
    
  my $vf_adaptor = Bio::EnsEMBL::DBSQL::MergedAdaptor->new(
    -species  => $self->adaptor()->db()->species, 
    -type     => "VariationFeature"
  );
  
  if( $vf_adaptor ) {
    return $vf_adaptor;
  }
  else {
    warning("Variation database must be attached to core database to " .
            "retrieve variation information" );
        
    return undef;
  }
}
*/
/*
=head2 get_all_VariationFeatures
    Args        : $so_terms [optional] - list of so_terms to limit the fetch to
    Description : Returns all germline variation features on this slice. This function will 
                  only work correctly if the variation database has been attached to the core 
                  database.
                  If $so_terms is specified, only variation features with a consequence type
                  that matches or is an ontological child of any of the supplied terms will
                  be returned
    ReturnType  : listref of Bio::EnsEMBL::Variation::VariationFeature
    Exceptions  : none
    Caller      : contigview, snpview
    Status      : Stable

=cut

sub get_all_VariationFeatures{
  my $self     = shift;

  if (my $vf_adaptor = $self->_get_VariationFeatureAdaptor) {
    return $vf_adaptor->fetch_all_by_Slice_SO_terms($self, @_);
  }
  else {
    return [];
  }
}

=head2 get_all_somatic_VariationFeatures

    Args        : $filter [optional]
    Description : Returns all somatic variation features on this slice. This function will only 
                  work correctly if the variation database has been attached to the core database.
    ReturnType  : listref of Bio::EnsEMBL::Variation::VariationFeature
    Exceptions  : none
    Status      : Stable

=cut

sub get_all_somatic_VariationFeatures {
  my $self = shift;
    
  if (my $vf_adaptor = $self->_get_VariationFeatureAdaptor) {
    return $vf_adaptor->fetch_all_somatic_by_Slice($self);
  }
  else{
    return [];
  }
}

=head2 get_all_somatic_VariationFeatures_by_source

    Args        : $source [optional]
    Description : Returns all somatic variation features, from a defined source name (e.g.'COSMIC'), 
                  on this slice. This function will only work correctly if the variation database
                  has been attached to the core database.
    ReturnType  : listref of Bio::EnsEMBL::Variation::VariationFeature
    Exceptions  : none
    Status      : Stable

=cut

sub get_all_somatic_VariationFeatures_by_source {
  my $self     = shift;
  my $source   = shift;
  my $constraint = (defined($source)) ? " s.name='$source' " : undef;
  
  if (my $vf_adaptor = $self->_get_VariationFeatureAdaptor) {
    return $vf_adaptor->fetch_all_somatic_by_Slice_constraint($self, $constraint);
  }
  else {
    return [];
  }
}


=head2 get_all_VariationFeatures_with_annotation

    Arg [1]     : $variation_feature_source [optional]
    Arg [2]     : $annotation_source [optional]
    Arg [3]     : $annotation_name [optional]
    Description : returns all germline variation features on this slice associated with a phenotype.
                  This function will only work correctly if the variation database has been
                  attached to the core database.
                  If $variation_feature_source is set only variations from that source
                  are retrieved.
                  If $annotation_source is set only variations whose annotations come from
                  $annotation_source will be retrieved.
                  If $annotation_name is set only variations with that annotation will be retrieved.
                  $annotation_name can be a phenotype's internal dbID.
    ReturnType  : listref of Bio::EnsEMBL::Variation::VariationFeature
    Exceptions  : none
    Caller      : contigview, snpview
    Status      : Stable

=cut

sub get_all_VariationFeatures_with_annotation{
  my $self = shift;
  my $source = shift;
  my $p_source = shift;
  my $annotation = shift;

  if (my $vf_adaptor = $self->_get_VariationFeatureAdaptor) {
    return $vf_adaptor->fetch_all_with_annotation_by_Slice($self, $source, $p_source, $annotation);
  }
  else {
    return [];
  }
}

=head2 get_all_somatic_VariationFeatures_with_annotation

    Arg [1]     : $variation_feature_source [optional]
    Arg [2]     : $annotation_source [optional]
    Arg [3]     : $annotation_name [optional]
    Description : returns all somatic variation features on this slice associated with a phenotype.
                  (see get_all_VariationFeatures_with_annotation for further documentation)
    ReturnType  : listref of Bio::EnsEMBL::Variation::VariationFeature
    Exceptions  : none
    Status      : Stable

=cut

sub get_all_somatic_VariationFeatures_with_annotation{
  my $self = shift;
  my $source = shift;
  my $p_source = shift;
  my $annotation = shift;

  if (my $vf_adaptor = $self->_get_VariationFeatureAdaptor) {
    return $vf_adaptor->fetch_all_somatic_with_annotation_by_Slice($self, $source, $p_source, $annotation);
  }
  else {
    return [] unless $vf_adaptor;
  }
}

=head2 get_all_VariationFeatures_by_VariationSet

    Arg [1]     : Bio::EnsEMBL:Variation::VariationSet $set
    Description :returns all variation features on this slice associated with a given set.
                 This function will only work correctly if the variation database has been
                 attached to the core database. 
    ReturnType : listref of Bio::EnsEMBL::Variation::VariationFeature
    Exceptions : none
    Caller     : contigview, snpview
    Status     : Stable

=cut

sub get_all_VariationFeatures_by_VariationSet {
  my $self = shift;
  my $set = shift;

  if (my $vf_adaptor = $self->_get_VariationFeatureAdaptor) {
    return $vf_adaptor->fetch_all_by_Slice_VariationSet($self, $set);  
  }
  else {
    return [];
  }
}

=head2 get_all_StructuralVariations

                Description: DEPRECATED. Use get_all_StructuralVariationFeatures instead

=cut

sub get_all_StructuralVariations{
  my $self = shift;
  my $source = shift;
        my $study = shift;
        my $sv_class = shift;
        
        deprecate('Use get_all_StructuralVariationFeatures() instead.');

        return $self->get_all_StructuralVariationFeatures($source,$sv_class);
}




=head2 get_all_CopyNumberVariantProbes

        Description: DEPRECATED. Use get_all_CopyNumberVariantProbeFeatures instead
        
=cut

sub get_all_CopyNumberVariantProbes {
        my $self = shift;
  my $source = shift;
        my $study = shift;
        
        deprecate('Use get_all_CopyNumberVariantProbeFeatures() instead.');
        
        return $self->get_all_CopyNumberVariantProbeFeatures($source);
}


sub _get_StructuralVariationFeatureAdaptor {
    
  my $self = shift;
    
  if(!$self->adaptor()) {
    warning('Cannot get structural variation features without attached adaptor');
    return undef;
  }
    
  my $svf_adaptor = Bio::EnsEMBL::DBSQL::MergedAdaptor->new(
    -species  => $self->adaptor()->db()->species, 
    -type     => "StructuralVariationFeature"
  );
  
  if( $svf_adaptor ) {
    return $svf_adaptor;
  }
  else {
    warning("Variation database must be attached to core database to " .
            "retrieve variation information" );
        
    return undef;
  }
}


=head2 get_all_StructuralVariationFeatures

    Arg[1]      : string $source [optional]
                Arg[2]      : int $include_evidence [optional]        
                Arg[3]      : string $sv_class (SO term) [optional]        
    Description : returns all structural variation features on this slice. This function will only work
                  correctly if the variation database has been attached to the core database.
                  If $source is set, only structural variation features with that source name will be 
                                                                        returned. By default, it only returns structural variant features which are not labelled 
                                                                        as "CNV_PROBE".
                                                                        If $include_evidence is set (i.e. $include_evidence=1), structural variation features from 
                                                            both structural variation (SV) and their supporting structural variations (SSV) will be 
                                                            returned. By default, it only returns features from structural variations (SV).
    ReturnType  : listref of Bio::EnsEMBL::Variation::StructuralVariationFeature
    Exceptions  : none
    Caller      : contigview, snpview, structural_variation_features
    Status      : Stable

=cut

sub get_all_StructuralVariationFeatures {
  my $self             = shift;
  my $source           = shift;
        my $include_evidence = shift;
        my $somatic          = shift;
        my $sv_class         = shift;
        my $constraint       = shift;
        
        my $operator = '';
        
  if (!defined($sv_class)) { 
                $sv_class = 'SO:0000051'; # CNV_PROBE
                $operator = '!'; # All but CNV_PROBE
        }
        
        $somatic = (!defined($somatic) || !$somatic) ? 0 : 1;
        
        my $svf_adaptor = $self->_get_StructuralVariationFeatureAdaptor;
        
        my $variation_db = $self->adaptor->db->get_db_adaptor('variation');
        
        # Get the attrib_id
        my $at_adaptor = $variation_db->get_AttributeAdaptor;
        my $SO_term   = $at_adaptor->SO_term_for_SO_accession($sv_class);
        my $attrib_id = $at_adaptor->attrib_id_for_type_value('SO_term',$SO_term);

        if (!$attrib_id) {
                warning("The Sequence Ontology accession number is not found in the database");
          return [];
        }
        
        # Get the structural variations features
  if( $svf_adaptor ) {
    $constraint .= qq{ AND } if ($constraint);
    $constraint .= qq{ svf.somatic=$somatic AND svf.class_attrib_id $operator=$attrib_id };
    $constraint .= qq{ AND svf.is_evidence=0 } if (!$include_evidence);

                if($source) {
      return $svf_adaptor->fetch_all_by_Slice_constraint($self, qq{$constraint AND s.name = '$source'});
    }else {
                        return $svf_adaptor->fetch_all_by_Slice_constraint($self, $constraint);
    }
  }
  else {
                warning("Variation database must be attached to core database to " .
                                                 "retrieve variation information" );
    return [];
  }
}

=head2 get_all_StructuralVariationFeatures_by_size_range
    Arg[1]      : int $size_min (minimum size of the structural variant)
    Arg[2]      : int $size_max (maximum size of the structural variant) [optional]
    Arg[1]      : int minimum size of the structural variant
    Arg[3]      : string $source [optional]
                Arg[4]      : int $include_evidence [optional]        
                Arg[5]      : string $sv_class (SO term) [optional]        
    Description : returns all structural variation features overlapping this slice with a size greater than the minimum size 
                              defined in the first argument, and (optional) lesser than the maximun size defined in the second argument. 
                  This function will only work correctly if the variation database has been attached to the core database.
                  If $source is set, only structural variation features with that source name will be 
                                                                        returned. By default, it only returns structural variant features which are not labelled 
                                                                        as "CNV_PROBE".
                                                                        If $include_evidence is set (i.e. $include_evidence=1), structural variation features from 
                                                            both structural variation (SV) and their supporting structural variations (SSV) will be 
                                                            returned. By default, it only returns features from structural variations (SV).
    ReturnType  : listref of Bio::EnsEMBL::Variation::StructuralVariationFeature
    Exceptions  : none
    Caller      : contigview, snpview, structural_variation_features
    Status      : Stable

=cut

sub get_all_StructuralVariationFeatures_by_size_range {
        my $self             = shift;
        my $size_min         = shift;
        my $size_max         = shift;
        my $source           = shift;
        my $include_evidence = shift;
        my $somatic          = shift;
        my $sv_class         = shift;
        
        my $constraint = qq{svf.seq_region_end-svf.seq_region_start>=$size_min};
           $constraint .= qq{ AND svf.seq_region_end-svf.seq_region_start<$size_max } if (defined $size_max);
        
        return $self->get_all_StructuralVariationFeatures($source,$include_evidence,$somatic,$sv_class,$constraint);
}

=head2 get_all_StructuralVariationFeatures_by_VariationSet

    Arg [1]     : Bio::EnsEMBL:Variation::VariationSet $set
    Description :returns all structural variation features on this slice associated with a 
                 given set.
                 This function will only work correctly if the variation database has been
                 attached to the core database. 
    ReturnType : listref of Bio::EnsEMBL::Variation::StructuralVariationFeature
    Exceptions : none
    Caller     : contigview, snpview
    Status     : Stable

=cut

sub get_all_StructuralVariationFeatures_by_VariationSet {
  my $self = shift;
  my $set = shift;

  if (my $svf_adaptor = $self->_get_StructuralVariationFeatureAdaptor) {
    return $svf_adaptor->fetch_all_by_Slice_VariationSet($self, $set);  
  }
  else {
    return [];
  }
}


=head2 get_all_somatic_StructuralVariationFeatures

                Arg[1]      : string $source [optional]
                Arg[2]      : int $include_evidence [optional]
                Arg[3]      : string $sv_class (SO term) [optional]                
    Description : returns all somatic structural variation features on this slice. This function will only work
                  correctly if the variation database has been attached to the core database.
                  If $source is set, only somatic structural variation features with that source name will be 
                                                                        returned. By default, it only returns somatic structural variant features which are not labelled 
                                                                        as "CNV_PROBE".
                                                                        If $include_evidence is set (i.e. $include_evidence=1), structural variation features from 
                                                            both structural variation (SV) and their supporting structural variations (SSV) will be 
                                                            returned. By default, it only returns features from structural variations (SV).
    ReturnType  : listref of Bio::EnsEMBL::Variation::StructuralVariationFeature
    Exceptions  : none
    Caller      : contigview, snpview, structural_variation_features
    Status      : Stable

=cut

sub get_all_somatic_StructuralVariationFeatures {
        my $self   = shift;
        my $source = shift;
  my $include_evidence = shift;
        my $sv_class = shift;
        
        return $self->get_all_StructuralVariationFeatures($source,$include_evidence,1,$sv_class);
}


=head2 get_all_CopyNumberVariantProbeFeatures

    Arg[1]      : string $source [optional]
    Description : returns all copy number variant probes on this slice. This function will only work
                  correctly if the variation database has been attached to the core database.
                  If $source is set, only CNV probes with that source name will be returned.
                  If $study is set, only CNV probes of that study will be returned.
    ReturnType  : listref of Bio::EnsEMBL::Variation::StructuralVariationFeature
    Exceptions  : none
    Caller      : contigview, snpview, structural_variation_feature
    Status      : At Risk

=cut

sub get_all_CopyNumberVariantProbeFeatures {
        my $self   = shift;
  my $source = shift;
        
        return $self->get_all_StructuralVariationFeatures($source,0,0,'SO:0000051');
}


=head2 get_all_VariationFeatures_by_Population

  Arg [1]    : Bio::EnsEMBL::Variation::Population
  Arg [2]         : $minimum_frequency (optional)
  Example    : $pop = $pop_adaptor->fetch_by_dbID(659);
               @vfs = @{$slice->get_all_VariationFeatures_by_Population(
                 $pop,$slice)};
  Description: Retrieves all variation features in a slice which are stored for
                           a specified population. If $minimum_frequency is supplied, only
                           variations with a minor allele frequency (MAF) greater than
                           $minimum_frequency will be returned.
  Returntype : listref of Bio::EnsEMBL::Variation::VariationFeature
  Exceptions : throw on incorrect argument
  Caller     : general
  Status     : At Risk

=cut

sub get_all_VariationFeatures_by_Population {
  my $self = shift;

  if (my $vf_adaptor = $self->_get_VariationFeatureAdaptor) {
    return $vf_adaptor->fetch_all_by_Slice_Population($self, @_);
  }
  else {
    return [];
  }
}




=head2 get_all_IndividualSlice

    Args        : none
    Example     : my $individualSlice = $slice->get_by_Population($population);
    Description : Gets the specific Slice for all the individuls in the population
    ReturnType  : listref of Bio::EnsEMB::IndividualSlice
    Exceptions  : none
    Caller      : general

=cut

sub get_all_IndividualSlice{
    my $self = shift;

    my $individualSliceFactory = Bio::EnsEMBL::IndividualSliceFactory->new(
                                                                           -START   => $self->{'start'},
                                                                           -END     => $self->{'end'},
                                                                           -STRAND  => $self->{'strand'},
                                                                           -ADAPTOR => $self->adaptor(),
                                                                           -SEQ_REGION_NAME => $self->{'seq_region_name'},
                                                                           -SEQ_REGION_LENGTH => $self->{'seq_region_length'},
                                                                           -COORD_SYSTEM    => $self->{'coord_system'},
                                                                           );
    return $individualSliceFactory->get_all_IndividualSlice();
}

=head2 get_by_Individual

    Arg[1]      : Bio::EnsEMBL::Variation::Individual $individual
    Example     : my $individualSlice = $slice->get_by_Individual($individual);
    Description : Gets the specific Slice for the individual
    ReturnType  : Bio::EnsEMB::IndividualSlice
    Exceptions  : none
    Caller      : general

=cut

sub get_by_Individual{
    my $self = shift;
    my $individual = shift;

    return Bio::EnsEMBL::IndividualSlice->new(
                                          -START   => $self->{'start'},
                                          -END     => $self->{'end'},
                                          -STRAND  => $self->{'strand'},
                                          -ADAPTOR => $self->adaptor(),
#                                          -SEQ     => $self->{'seq'},
                                          -SEQ_REGION_NAME => $self->{'seq_region_name'},
                                          -SEQ_REGION_LENGTH => $self->{'seq_region_length'},
                                          -COORD_SYSTEM    => $self->{'coord_system'},
                                          -INDIVIDUAL     => $individual);

}



=head2 get_by_strain

    Arg[1]      : string $strain
    Example     : my $strainSlice = $slice->get_by_strain($strain);
    Description : Gets the specific Slice for the strain
    ReturnType  : Bio::EnsEMB::StrainSlice
    Exceptions  : none
    Caller      : general

=cut

sub get_by_strain{
    my $self = shift;
    my $strain_name = shift;

    return Bio::EnsEMBL::StrainSlice->new(
                                          -START   => $self->{'start'},
                                          -END     => $self->{'end'},
                                          -STRAND  => $self->{'strand'},
                                          -ADAPTOR => $self->adaptor(),
                                          -SEQ     => $self->{'seq'},
                                          -SEQ_REGION_NAME => $self->{'seq_region_name'},
                                          -SEQ_REGION_LENGTH => $self->{'seq_region_length'},
                                          -COORD_SYSTEM    => $self->{'coord_system'},
                                          -STRAIN_NAME     => $strain_name);

}

sub calculate_theta{
    my $self = shift;
    my $strains = shift;
    my $feature = shift; #optional parameter. Name of the feature in the Slice you want to calculate

    if(!$self->adaptor()) {
        warning('Cannot get variation features without attached adaptor');
        return 0;
    }
    my $variation_db = $self->adaptor->db->get_db_adaptor('variation');

    unless($variation_db) {
        warning("Variation database must be attached to core database to " .
                "retrieve variation information" );
        return 0;
    }

    #need to get coverage regions for the slice in the different strains
    my $coverage_adaptor = $variation_db->get_ReadCoverageAdaptor;
    my $strain;
    my $differences = [];
    my $slices = [];
    if ($coverage_adaptor){
        my $num_strains = scalar(@{$strains}) +1;
        if (!defined $feature){
            #we want to calculate for the whole slice
            push @{$slices}, $self; #add the slice as the slice to calculate the theta value
        }
        else{
            #we have features, get the slices for the different features
            my $features = $self->get_all_Exons();
            map {push @{$slices},$_->feature_Slice} @{$features}; #add the slices of the features
        }
        my $length_regions = 0;
        my $snps = 0;
        my $theta = 0;
        my $last_position = 0;
        #get all the differences in the slice coordinates
        foreach my $strain_name (@{$strains}){
            my $strain = $self->get_by_strain($strain_name); #get the strainSlice for the strain

            my $results = $strain->get_all_differences_Slice;
            push @{$differences}, @{$results} if (defined $results);
        }
        #when we finish, we have, in max_level, the regions covered by all the sample
        #sort the differences by the genomic position
        my @differences_sorted = sort {$a->start <=> $b->start} @{$differences};
        foreach my $slice (@{$slices}){
            my $regions_covered = $coverage_adaptor->fetch_all_regions_covered($slice,$strains);
            if (defined $regions_covered){
                foreach my $range (@{$regions_covered}){
                    $length_regions += ($range->[1] - $range->[0]) + 1; #add the length of the genomic region
                    for (my $i = $last_position;$i<@differences_sorted;$i++){
                        if ($differences_sorted[$i]->start >= $range->[0] && $differences_sorted[$i]->end <= $range->[1]){
                            $snps++; #count differences in the region
                        }
                        elsif ($differences_sorted[$i]->end > $range->[1]){
                            $last_position = $i;
                            last;
                        }
                    }
                }
                #when all the ranges have been iterated, calculate rho
                #this is an intermediate variable called a in the formula
                #  a = sum i=2..strains 1/i-1
            }
        }
        my $a = _calculate_a($num_strains);
        $theta = $snps / ($a * $length_regions);
        return $theta;
    }
    else{
        return 0;
    }
}




sub _calculate_a{
    my $max_level = shift;

    my $a = 0;
    for (my $i = 2; $i <= $max_level+1;$i++){
        $a += 1/($i-1);
    }
    return $a;
}

sub calculate_pi{
    my $self = shift;
    my $strains = shift;
    my $feature = shift;

    if(!$self->adaptor()) {
        warning('Cannot get variation features without attached adaptor');
        return 0;
    }
    my $variation_db = $self->adaptor->db->get_db_adaptor('variation');

    unless($variation_db) {
        warning("Variation database must be attached to core database to " .
                "retrieve variation information" );
        return 0;
    }

    #need to get coverage regions for the slice in the different strains
    my $coverage_adaptor = $variation_db->get_ReadCoverageAdaptor;
    my $differences = [];
    my $slices = [];
    if ($coverage_adaptor){
        my $num_strains = scalar(@{$strains}) +1;
        if (!defined $feature){
            #we want to calculate for the whole slice
            push @{$slices}, $self; #add the slice as the slice to calculate the theta value
        }
        else{
            #we have features, get the slices for the different features
            my $features = $self->get_all_Exons();
            map {push @{$slices},$_->feature_Slice} @{$features}; #add the slices of the features
        }
        my @range_differences = ();
        my $pi = 0;
        my $regions = 0;
        my $last_position = 0; #last position visited in the sorted list of differences
        my $triallelic = 0;
        my $is_triallelic = 0;
        foreach my $slice (@{$slices}){
            foreach my $strain_name (@{$strains}){
                my $strain = $slice->get_by_strain($strain_name); #get the strainSlice for the strain
                my $results = $strain->get_all_differences_Slice;
                push @{$differences}, @{$results} if (defined $results);
            }
            my @differences_sorted = sort {$a->start <=> $b->start} @{$differences};

            my $regions_covered = $coverage_adaptor->fetch_all_regions_covered($slice,$strains);
            #when we finish, we have, in max_level, the regions covered by all the sample
            #sort the differences
            if (defined $regions_covered){
                foreach my $range (@{$regions_covered}){
                    for (my $i = $last_position;$i<@differences_sorted;$i++){
                        if ($differences_sorted[$i]->start >= $range->[0] && $differences_sorted[$i]->end <= $range->[1]){
                            #check wether it is the same region or different
                            if (!defined $range_differences[0] || ($differences_sorted[$i]->start == $range_differences[0]->start)){
                                if (defined $range_differences[0] && ($differences_sorted[$i]->allele_string ne $range_differences[0]->allele_string)){
                                    $is_triallelic = 1;
                                }
                                push @range_differences, $differences_sorted[$i];
                            }
                            else{
                                #new site, calc pi for the previous one
                                $pi += 2 * (@range_differences/($num_strains)) * ( 1 - (@range_differences/$num_strains));
                                if ($is_triallelic) {
                                    $triallelic++;
                                    $is_triallelic = 0;
                                }
                                $regions++;
                                @range_differences = ();
                                #and start a new range
                                push @range_differences, $differences_sorted[$i];
                            }
                        }
                        elsif ($differences_sorted[$i]->end > $range->[1]){
                            $last_position = $i;
                            last;
                        }
                    }
                    #calculate pi for last site, if any
                    if (defined $range_differences[0]){
                        $pi += 2 * (@range_differences/$num_strains) * ( 1 - (@range_differences/$num_strains));
                        $regions++;
                    }
                }
            }
            $pi = $pi / $regions; #calculate average pi
            print "Regions with variations in region $regions and triallelic $triallelic\n\n";
        }
        return $pi;
    }
    else{
        return 0;
    }

}





=head2 get_all_genotyped_VariationFeatures

    Args       : none
    Function   : returns all variation features on this slice that have been genotyped. This function will only work
                correctly if the variation database has been attached to the core database.
    ReturnType : listref of Bio::EnsEMBL::Variation::VariationFeature
    Exceptions : none
    Caller     : contigview, snpview, ldview
    Status     : At Risk
               : Variation database is under development.

=cut

sub get_all_genotyped_VariationFeatures{
  my $self = shift;

  if( my $vf_adaptor = $self->_get_VariationFeatureAdaptor) {
    return $vf_adaptor->fetch_all_genotyped_by_Slice($self);
  } 
  else {
    return [];
  }
}


=head2 get_all_SNPs

 Description: DEPRECATED. Use get_all_VariationFeatures insted

=cut

sub get_all_SNPs {
  my $self = shift;

  deprecate('Use get_all_VariationFeatures() instead.');

  my $snps;
  my $vf = $self->get_all_genotyped_VariationFeatures();
  if( $vf->[0] ) {
      #necessary to convert the VariationFeatures into SNP objects
      foreach my $variation_feature (@{$vf}){
          push @{$snps},$variation_feature->convert_to_SNP();
      }
      return $snps;
  } else {
    return [];
  }
}

=head2 get_all_genotyped_SNPs

  Description   : DEPRECATED. Use get_all_genotyped_VariationFeatures insted

=cut

sub get_all_genotyped_SNPs {
  my $self = shift;

  deprecate("Use get_all_genotyped_VariationFeatures instead");
  my $vf = $self->get_all_genotyped_VariationFeatures;
  my $snps;
  if ($vf->[0]){
      foreach my $variation_feature (@{$vf}){
          push @{$snps},$variation_feature->convert_to_SNP();
      }
      return $snps;
  } else {
      return [];
  }
}

sub get_all_SNPs_transcripts {
  my $self = shift;

  deprecate("DEPRECATED");

  return [];

}
*/



/*
=head2 get_all_Genes

  Arg [1]    : (optional) string $logic_name
               The name of the analysis used to generate the genes to retrieve
  Arg [2]    : (optional) string $dbtype
               The dbtype of genes to obtain.  This assumes that the db has
               been added to the DBAdaptor under this name (using the
               DBConnection::add_db_adaptor method).
  Arg [3]    : (optional) boolean $load_transcripts
               If set to true, transcripts will be loaded immediately rather
               than being lazy-loaded on request.  This will result in a
               significant speed up if the Transcripts and Exons are going to
               be used (but a slow down if they are not).
  Arg [4]    : (optional) string $source
               The source of the genes to retrieve.
  Arg [5]    : (optional) string $biotype
               The biotype of the genes to retrieve.
  Example    : @genes = @{$slice->get_all_Genes};
  Description: Retrieves all genes that overlap this slice.
  Returntype : listref of Bio::EnsEMBL::Genes
  Exceptions : none
  Caller     : none
  Status     : Stable

=cut
*/

Vector *Slice_getAllGenes(Slice *slice, char *logicName, char *dbType, int loadTranscripts, char *source, char *bioType) {
  if (!Slice_getAdaptor(slice)) {
    fprintf(stderr, "Warning: Cannot get Genes without attached adaptor");
    return Vector_new();
  }

// Note: Gene did this db type selection in an even more peculiar way using the registry - I'm just going to use the same way as other feature types
  DBAdaptor *db = Slice_getSelectedDBAdaptor(slice, dbType);

  if (!db) {
    fprintf(stderr, "Warning: Don't have db %s returning empty list\n", dbType);
    return Vector_new();
  }

  GeneAdaptor *ga = DBAdaptor_getGeneAdaptor(db);

  return GeneAdaptor_fetchAllBySlice(ga, slice, logicName, loadTranscripts, source, bioType);
}

/*
=head2 get_all_Genes_by_type

  Arg [1]    : string $type
               The biotype of genes wanted.
  Arg [2]    : (optional) string $logic_name
  Arg [3]    : (optional) boolean $load_transcripts
               If set to true, transcripts will be loaded immediately rather
               than being lazy-loaded on request.  This will result in a
               significant speed up if the Transcripts and Exons are going to
               be used (but a slow down if they are not).
  Example    : @genes = @{$slice->get_all_Genes_by_type('protein_coding',
               'ensembl')};
  Description: Retrieves genes that overlap this slice of biotype $type.
               This is primarily used by the genebuilding code when several
               biotypes of genes are used.

               The logic name is the analysis of the genes that are retrieved.
               If not provided all genes will be retrieved instead.

  Returntype : listref of Bio::EnsEMBL::Genes
  Exceptions : none
  Caller     : genebuilder, general
  Status     : Stable

=cut
*/

/*
Vector *Slice_getAllGenesByType(Slice *slice, char *type, char *logicName, int loadTranscripts) {
  if (!Slice_getAdaptor(slice)) {
    fprintf(stderr, "Warning: Cannot get Genes without attached adaptor");
    return Vector_new();
  }

  return Slice_getAllGenes(slice, logicName, NULL, loadTranscripts, NULL, type);
}
*/


/*
=head2 get_all_Genes_by_source

  Arg [1]    : string source
  Arg [2]    : (optional) boolean $load_transcripts
               If set to true, transcripts will be loaded immediately rather
               than being lazy-loaded on request.  This will result in a
               significant speed up if the Transcripts and Exons are going to
               be used (but a slow down if they are not).
  Example    : @genes = @{$slice->get_all_Genes_by_source('ensembl')};
  Description: Retrieves genes that overlap this slice of source $source.

  Returntype : listref of Bio::EnsEMBL::Genes
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/

/*
Vector *Slice_getAllGenesBySource(Slice *slice, char *source, int loadTranscripts) {
  if (!Slice_getAdaptor(slice)) {
    fprintf(stderr, "Warning: Cannot get Genes without attached adaptor");
    return Vector_new();
  }

  return Slice_getAllGenes(slice, NULL, NULL, loadTranscripts, source, NULL);
}
*/

/*
=head2 get_all_Transcripts

  Arg [1]    : (optional) boolean $load_exons
               If set to true exons will not be lazy-loaded but will instead
               be loaded right away.  This is faster if the exons are
               actually going to be used right away.
  Arg [2]    : (optional) string $logic_name
               the logic name of the type of features to obtain
  Arg [3]    : (optional) string $db_type
  Example    : @transcripts = @{$slice->get_all_Transcripts)_};
  Description: Gets all transcripts which overlap this slice.  If you want to
               specify a particular analysis or type, then you are better off
               using get_all_Genes or get_all_Genes_by_type and iterating
               through the transcripts of each gene.
  Returntype : reference to a list of Bio::EnsEMBL::Transcripts
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/

/*
Vector *Slice_getAllTranscripts(Slice *slice, int loadExons, char *logicName, char *dbType) {
  if (!Slice_getAdaptor(slice)) {
    fprintf(stderr, "Warning: Cannot get Transcripts without attached adaptor");
    return Vector_new();
  }

// Note: Transcript did this db type selection in an even more peculiar way using the registry - I'm just going to use the same way as other feature types
  DBAdaptor *db = Slice_getSelectedDBAdaptor(slice, dbType);

  if (!db) {
    fprintf(stderr, "Warning: Don't have db %s returning empty list\n", dbType);
    return Vector_new();
  }

  TranscriptAdaptor *ta = DBAdaptor_getTranscriptAdaptor(db);

  return TranscriptAdaptor_fetchAllBySlice(ta, slice, loadExons, logicName);
}
*/


/*
=head2 get_all_Exons

  Arg [1]    : none
  Example    : @exons = @{$slice->get_all_Exons};
  Description: Gets all exons which overlap this slice.  Note that these exons
               will not be associated with any transcripts, so this may not
               be terribly useful.
  Returntype : reference to a list of Bio::EnsEMBL::Exons
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/

/*
Vector *Slice_getAllExons(Slice *slice, char *dbType) {
  if (!Slice_getAdaptor(slice)) {
    fprintf(stderr, "Warning: Cannot get Exons without attached adaptor");
    return Vector_new();
  }

// Note: Exon didn't do dbType - added it to be consistent
  DBAdaptor *db = Slice_getSelectedDBAdaptor(slice, dbType);

  if (!db) {
    fprintf(stderr, "Warning: Don't have db %s returning empty list\n", dbType);
    return Vector_new();
  }

  ExonAdaptor *ea = DBAdaptor_getExonAdaptor(db);

  return ExonAdaptor_fetchAllBySlice(ea, slice);
}
*/


/*
=head2 get_all_QtlFeatures

  Args       : none
  Example    : none
  Description: returns overlapping QtlFeatures
  Returntype : listref Bio::EnsEMBL::Map::QtlFeature
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_QtlFeatures {
  my $self = shift;

  if(!$self->adaptor()) {
    warning('Cannot get QtlFeatures without attached adaptor');
    return [];
  }

  my $qfAdaptor;
  if( $self->adaptor()) {
    $qfAdaptor = $self->adaptor()->db()->get_QtlFeatureAdaptor();
  } else {
    return [];
  }

  return $qfAdaptor->fetch_all_by_Slice_constraint( $self );
}




=head2 get_all_KaryotypeBands

  Arg [1]    : none
  Example    : @kary_bands = @{$slice->get_all_KaryotypeBands};
  Description: Retrieves the karyotype bands which this slice overlaps.
  Returntype : listref oif Bio::EnsEMBL::KaryotypeBands
  Exceptions : none
  Caller     : general, contigview
  Status     : Stable

=cut

sub get_all_KaryotypeBands {
  my ($self) = @_;

  if(!$self->adaptor()) {
    warning('Cannot get KaryotypeBands without attached adaptor');
    return [];
  }

  my $kadp = $self->adaptor->db->get_KaryotypeBandAdaptor();
  return $kadp->fetch_all_by_Slice($self);
}

*/


/*
=head2 get_repeatmasked_seq

  Arg [1]    : listref of strings $logic_names (optional)
  Arg [2]    : int $soft_masking_enable (optional)
  Arg [3]    : hash reference $not_default_masking_cases (optional, default is {})
               The values are 0 or 1 for hard and soft masking respectively
               The keys of the hash should be of 2 forms
               "repeat_class_" . $repeat_consensus->repeat_class,
                e.g. "repeat_class_SINE/MIR"
               "repeat_name_" . $repeat_consensus->name
                e.g. "repeat_name_MIR"
               depending on which base you want to apply the not default
               masking either the repeat_class or repeat_name. Both can be
               specified in the same hash at the same time, but in that case,
               repeat_name setting has priority over repeat_class. For example,
               you may have hard masking as default, and you may want soft
               masking of all repeat_class SINE/MIR, but repeat_name AluSp
               (which are also from repeat_class SINE/MIR).
               Your hash will be something like {"repeat_class_SINE/MIR" => 1,
                                                 "repeat_name_AluSp" => 0}
  Example    : $rm_slice = $slice->get_repeatmasked_seq();
               $softrm_slice = $slice->get_repeatmasked_seq(['RepeatMask'],1);
  Description: Returns Bio::EnsEMBL::Slice that can be used to create repeat
               masked sequence instead of the regular sequence.
               Sequence returned by this new slice will have repeat regions
               hardmasked by default (sequence replaced by N) or
               or soft-masked when arg[2] = 1 (sequence in lowercase)
               Will only work with database connection to get repeat features.
  Returntype : Bio::EnsEMBL::RepeatMaskedSlice
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/

/* NIY: Should do this 
sub get_repeatmasked_seq {
    my ($self,$logic_names,$soft_mask,$not_default_masking_cases) = @_;

    return Bio::EnsEMBL::RepeatMaskedSlice->new
      (-START   => $self->{'start'},
       -END     => $self->{'end'},
       -STRAND  => $self->{'strand'},
       -ADAPTOR => $self->adaptor(),
       -SEQ     => $self->{'seq'},
       -SEQ_REGION_NAME => $self->{'seq_region_name'},
       -SEQ_REGION_LENGTH => $self->{'seq_region_length'},
       -COORD_SYSTEM    => $self->{'coord_system'},
       -REPEAT_MASK     => $logic_names,
       -SOFT_MASK       => $soft_mask,
       -NOT_DEFAULT_MASKING_CASES => $not_default_masking_cases);
}
*/


/*

=head2 _mask_features

  Arg [1]    : reference to a string $dnaref
  Arg [2]    : array_ref $repeats
               reference to a list Bio::EnsEMBL::RepeatFeature
               give the list of coordinates to replace with N or with
               lower case
  Arg [3]    : int $soft_masking_enable (optional)
  Arg [4]    : hash reference $not_default_masking_cases (optional, default is {})
               The values are 0 or 1 for hard and soft masking respectively
               The keys of the hash should be of 2 forms
               "repeat_class_" . $repeat_consensus->repeat_class,
                e.g. "repeat_class_SINE/MIR"
               "repeat_name_" . $repeat_consensus->name
                e.g. "repeat_name_MIR"
               depending on which base you want to apply the not default masking either
               the repeat_class or repeat_name. Both can be specified in the same hash
               at the same time, but in that case, repeat_name setting has priority over
               repeat_class. For example, you may have hard masking as default, and
               you may want soft masking of all repeat_class SINE/MIR,
               but repeat_name AluSp (which are also from repeat_class SINE/MIR).
               Your hash will be something like {"repeat_class_SINE/MIR" => 1,
                                                 "repeat_name_AluSp" => 0}
  Example    : none
  Description: replaces string positions described in the RepeatFeatures
               with Ns (default setting), or with the lower case equivalent
               (soft masking).  The reference to a dna string which is passed
               is changed in place.
  Returntype : none
  Exceptions : none
  Caller     : seq
  Status     : Stable

=cut
*/

/* NIY: Should do this 
sub _mask_features {
  my ($self,$dnaref,$repeats,$soft_mask,$not_default_masking_cases) = @_;

  $soft_mask = 0 unless (defined $soft_mask);
  $not_default_masking_cases = {} unless (defined $not_default_masking_cases);

  // explicit CORE::length call, to avoid any confusion with the Slice
  // length method
  my $dnalen = CORE::length($$dnaref);

 REP:foreach my $old_f (@{$repeats}) {
    my $f = $old_f->transfer( $self );
    my $start  = $f->start;
    my $end    = $f->end;
    my $length = ($end - $start) + 1;

    // check if we get repeat completely outside of expected slice range
    if ($end < 1 || $start > $dnalen) {
      //# warning("Unexpected: Repeat completely outside slice coordinates.");
      next REP;
    }

    // repeat partly outside slice range, so correct
    // the repeat start and length to the slice size if needed
    if ($start < 1) {
      $start = 1;
      $length = ($end - $start) + 1;
    }

    # repeat partly outside slice range, so correct
    # the repeat end and length to the slice size if needed
    if ($end > $dnalen) {
      $end = $dnalen;
      $length = ($end - $start) + 1;
    }

    $start--;

    my $padstr;
    // if we decide to define masking on the base of the repeat_type, we'll need
    // to add the following, and the other commented line few lines below.
    //# my $rc_type = "repeat_type_" . $f->repeat_consensus->repeat_type;
    my $rc_class = "repeat_class_" . $f->repeat_consensus->repeat_class;
    my $rc_name = "repeat_name_" . $f->repeat_consensus->name;

    my $masking_type;
    //# $masking_type = $not_default_masking_cases->{$rc_type} if (defined $not_default_masking_cases->{$rc_type});
    $masking_type = $not_default_masking_cases->{$rc_class} if (defined $not_default_masking_cases->{$rc_class});
    $masking_type = $not_default_masking_cases->{$rc_name} if (defined $not_default_masking_cases->{$rc_name});

    $masking_type = $soft_mask unless (defined $masking_type);

    if ($masking_type) {
      $padstr = lc substr ($$dnaref,$start,$length);
    } else {
      $padstr = 'N' x $length;
    }
    substr ($$dnaref,$start,$length) = $padstr;
  }
}
*/

/*

=head2 get_all_SearchFeatures

  Arg [1]    : scalar $ticket_ids
  Example    : $slice->get_all_SearchFeatures('BLA_KpUwwWi5gY');
  Description: Retreives all search features for stored blast
               results for the ticket that overlap this slice
  Returntype : listref of Bio::EnsEMBL::SeqFeatures
  Exceptions : none
  Caller     : general (webby!)
  Status     : Stable

=cut

sub get_all_SearchFeatures {
  my $self = shift;
  my $ticket = shift;
  local $_;
  unless($ticket) {
    throw("ticket argument is required");
  }

  if(!$self->adaptor()) {
    warning("Cannot get SearchFeatures without an attached adaptor");
    return [];
  }

  my $sfa = $self->adaptor()->db()->get_db_adaptor('blast');

  my $offset = $self->start-1;

  my $features = $sfa ? $sfa->get_all_SearchFeatures($ticket, $self->seq_region_name, $self->start, $self->end) : [];

  foreach( @$features ) {
    $_->start( $_->start - $offset );
    $_->end(   $_->end   - $offset );
  };
  return $features;

}
*/

/*
=head2 get_all_AssemblyExceptionFeatures

  Arg [1]    : string $set (optional)
  Example    : $slice->get_all_AssemblyExceptionFeatures();
  Description: Retreives all misc features which overlap this slice. If
               a set code is provided only features which are members of
               the requested set are returned.
  Returntype : listref of Bio::EnsEMBL::AssemblyExceptionFeatures
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
sub get_all_AssemblyExceptionFeatures {
  my $self = shift;
  my $misc_set = shift;

  my $adaptor = $self->adaptor();

  if(!$adaptor) {
    warning('Cannot retrieve features without attached adaptor.');
    return [];
  }

  my $aefa = $adaptor->db->get_AssemblyExceptionFeatureAdaptor();

  return $aefa->fetch_all_by_Slice($self);
}
*/

/*

=head2 get_all_MiscFeatures

  Arg [1]    : string $set (optional)
  Arg [2]    : string $database (optional)
  Example    : $slice->get_all_MiscFeatures('cloneset');
  Description: Retreives all misc features which overlap this slice. If
               a set code is provided only features which are members of
               the requested set are returned.
  Returntype : listref of Bio::EnsEMBL::MiscFeatures
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_MiscFeatures {
  my $self = shift;
  my $misc_set = shift;
  my $dbtype = shift;
  my $msa;

  my $adaptor = $self->adaptor();
  if(!$adaptor) {
    warning('Cannot retrieve features without attached adaptor.');
    return [];
  }

  my $mfa;
  if($dbtype) {
    my $db = $reg->get_db($adaptor->db(), $dbtype);
    if(defined($db)){
      $mfa = $reg->get_adaptor( lc($db->species()), $db->group(), "miscfeature" );
    } else{
      $mfa = $reg->get_adaptor( $adaptor->db()->species(), $dbtype, "miscfeature" );
    }
    if(!defined $mfa) {
      warning( "$dbtype misc features not available" );
      return [];
    }
  } else {
    $mfa =  $adaptor->db->get_MiscFeatureAdaptor();
  }

  if($misc_set) {
    return $mfa->fetch_all_by_Slice_and_set_code($self,$misc_set);
  }

  return $mfa->fetch_all_by_Slice($self);
}

=head2 get_all_MarkerFeatures

  Arg [1]    : (optional) string logic_name
               The logic name of the marker features to retrieve
  Arg [2]    : (optional) int $priority
               Lower (exclusive) priority bound of the markers to retrieve
  Arg [3]    : (optional) int $map_weight
               Upper (exclusive) priority bound of the markers to retrieve
  Example    : my @markers = @{$slice->get_all_MarkerFeatures(undef,50, 2)};
  Description: Retrieves all markers which lie on this slice fulfilling the
               specified map_weight and priority parameters (if supplied).
  Returntype : reference to a list of Bio::EnsEMBL::MarkerFeatures
  Exceptions : none
  Caller     : contigview, general
  Status     : Stable

=cut

sub get_all_MarkerFeatures {
  my ($self, $logic_name, $priority, $map_weight) = @_;

  if(!$self->adaptor()) {
    warning('Cannot retrieve MarkerFeatures without attached adaptor.');
    return [];
  }

  my $ma = $self->adaptor->db->get_MarkerFeatureAdaptor;

  my $feats = $ma->fetch_all_by_Slice_and_priority($self,
                                              $priority,
                                              $map_weight,
                                              $logic_name);
  return $feats;
}


=head2 get_MarkerFeatures_by_Name

  Arg [1]    : string marker Name
               The name (synonym) of the marker feature(s) to retrieve
  Example    : my @markers = @{$slice->get_MarkerFeatures_by_Name('z1705')};
  Description: Retrieves all markers with this ID
  Returntype : reference to a list of Bio::EnsEMBL::MarkerFeatures
  Exceptions : none
  Caller     : contigview, general
  Status     : Stable

=cut

sub get_MarkerFeatures_by_Name {
  my ($self, $name) = @_;

  if(!$self->adaptor()) {
    warning('Cannot retrieve MarkerFeatures without attached adaptor.');
    return [];
  }

  my $ma = $self->adaptor->db->get_MarkerFeatureAdaptor;

  my $feats = $ma->fetch_all_by_Slice_and_MarkerName($self, $name);
  return $feats;
}


=head2 get_all_compara_DnaAlignFeatures

  Arg [1]    : string $qy_species
               The name of the species to retrieve similarity features from
  Arg [2]    : string $qy_assembly
               The name of the assembly to retrieve similarity features from
  Arg [3]    : string $type
               The type of the alignment to retrieve similarity features from
  Arg [4]    : <optional> compara dbadptor to use.
  Example    : $fs = $slc->get_all_compara_DnaAlignFeatures('Mus musculus',
                                                            'MGSC3',
                                                            'WGA');
  Description: Retrieves a list of DNA-DNA Alignments to the species specified
               by the $qy_species argument.
               The compara database must be attached to the core database
               for this call to work correctly.  As well the compara database
               must have the core dbadaptors for both this species, and the
               query species added to function correctly.
  Returntype : reference to a list of Bio::EnsEMBL::DnaDnaAlignFeatures
  Exceptions : warning if compara database is not available
  Caller     : contigview
  Status     : Stable

=cut

sub get_all_compara_DnaAlignFeatures {
  my ($self, $qy_species, $qy_assembly, $alignment_type, $compara_db) = @_;

  if(!$self->adaptor()) {
    warning("Cannot retrieve DnaAlignFeatures without attached adaptor");
    return [];
  }

  unless($qy_species && $alignment_type # && $qy_assembly
  ) {
    throw("Query species and assembly and alignmemt type arguments are required");
  }

  if(!defined($compara_db)){
    $compara_db = Bio::EnsEMBL::Registry->get_DBAdaptor("compara", "compara");
  }
  unless($compara_db) {
    warning("Compara database must be attached to core database or passed ".
            "as an argument to " .
            "retrieve compara information");
    return [];
  }

  my $dafa = $compara_db->get_DnaAlignFeatureAdaptor;
  return $dafa->fetch_all_by_Slice($self, $qy_species, $qy_assembly, $alignment_type);
}

=head2 get_all_compara_Syntenies

  Arg [1]    : string $query_species e.g. "Mus_musculus" or "Mus musculus"
  Arg [2]    : string $method_link_type, default is "SYNTENY"
  Arg [3]    : <optional> compara dbadaptor to use.
  Description: gets all the compara syntenyies for a specfic species
  Returns    : arrayref of Bio::EnsEMBL::Compara::SyntenyRegion
  Status     : Stable

=cut

sub get_all_compara_Syntenies {
  my ($self, $qy_species, $method_link_type, $compara_db) = @_;

  if(!$self->adaptor()) {
    warning("Cannot retrieve features without attached adaptor");
    return [];
  }

  unless($qy_species) {
    throw("Query species and assembly arguments are required");
  }

  unless (defined $method_link_type) {
    $method_link_type = "SYNTENY";
  }

  if(!defined($compara_db)){
    $compara_db = Bio::EnsEMBL::Registry->get_DBAdaptor("compara", "compara");
  }
  unless($compara_db) {
    warning("Compara database must be attached to core database or passed ".
            "as an argument to " .
            "retrieve compara information");
    return [];
  }
  my $gdba = $compara_db->get_GenomeDBAdaptor();
  my $mlssa = $compara_db->get_MethodLinkSpeciesSetAdaptor();
  my $dfa = $compara_db->get_DnaFragAdaptor();
  my $sra = $compara_db->get_SyntenyRegionAdaptor();

  my $this_gdb = $gdba->fetch_by_core_DBAdaptor($self->adaptor()->db());
  my $query_gdb = $gdba->fetch_by_registry_name($qy_species);
  my $mlss = $mlssa->fetch_by_method_link_type_GenomeDBs($method_link_type, [$this_gdb, $query_gdb]);

  my $cs = $self->coord_system()->name();
  my $sr = $self->seq_region_name();
  my ($dnafrag) = @{$dfa->fetch_all_by_GenomeDB_region($this_gdb, $cs, $sr)};
  return $sra->fetch_all_by_MethodLinkSpeciesSet_DnaFrag($mlss, $dnafrag, $self->start, $self->end);
}

=head2 get_all_Haplotypes

  Arg [1]    : (optional) boolean $lite_flag
               if true lightweight haplotype objects are used
  Example    : @haplotypes = $slice->get_all_Haplotypes;
  Description: Retrieves all of the haplotypes on this slice.  Only works
               if the haplotype adaptor has been attached to the core adaptor
               via $dba->add_db_adaptor('haplotype', $hdba);
  Returntype : listref of Bio::EnsEMBL::External::Haplotype::Haplotypes
  Exceptions : warning is Haplotype database is not available
  Caller     : contigview, general
  Status     : Stable

=cut

sub get_all_Haplotypes {
  my($self, $lite_flag) = @_;

  if(!$self->adaptor()) {
    warning("Cannot retrieve features without attached adaptor");
    return [];
  }

  my $haplo_db = $self->adaptor->db->get_db_adaptor('haplotype');

  unless($haplo_db) {
    warning("Haplotype database must be attached to core database to " .
                "retrieve haplotype information" );
    return [];
  }

  my $haplo_adaptor = $haplo_db->get_HaplotypeAdaptor;

  my $haplotypes = $haplo_adaptor->fetch_all_by_Slice($self, $lite_flag);

  return $haplotypes;
}


sub get_all_DASFactories {
   my $self = shift;
   return [ $self->adaptor()->db()->_each_DASFeatureFactory ];
}

sub get_all_DASFeatures_dsn {
   my ($self, $source_type, $dsn) = @_;

  if(!$self->adaptor()) {
    warning("Cannot retrieve features without attached adaptor");
    return [];
  }
  my @X = grep { $_->adaptor->dsn eq $dsn } $self->adaptor()->db()->_each_DASFeatureFactory;

  return [ $X[0]->fetch_all_Features( $self, $source_type ) ];
}

=head2 get_all_DAS_Features

  Arg [1]    : none
  Example    : $features = $slice->get_all_DASFeatures;
  Description: Retreives a hash reference to a hash of DAS feature
               sets, keyed by the DNS, NOTE the values of this hash
               are an anonymous array containing:
                (1) a pointer to an array of features;
                (2) a pointer to the DAS stylesheet
  Returntype : hashref of Bio::SeqFeatures
  Exceptions : ?
  Caller     : webcode
  Status     : Stable

=cut
sub get_all_DAS_Features{
  my ($self) = @_;

  $self->{_das_features} ||= {}; # Cache
  $self->{_das_styles} ||= {}; # Cache
  $self->{_das_segments} ||= {}; # Cache
  my %das_features;
  my %das_styles;
  my %das_segments;
  my $slice = $self;

  foreach my $dasfact( @{$self->get_all_DASFactories} ){
    my $dsn = $dasfact->adaptor->dsn;
    my $name = $dasfact->adaptor->name;
#    my $type = $dasfact->adaptor->type;
    my $url = $dasfact->adaptor->url;

 my ($type) = $dasfact->adaptor->mapping;
 if (ref $type eq 'ARRAY') {
   $type = shift @$type;
 }
 $type ||= $dasfact->adaptor->type;
    # Construct a cache key : SOURCE_URL/TYPE
    # Need the type to handle sources that serve multiple types of features

    my $key = join('/', $name, $type);
    if( $self->{_das_features}->{$key} ){ # Use cached
        $das_features{$name} = $self->{_das_features}->{$key};
        $das_styles{$name} = $self->{_das_styles}->{$key};
        $das_segments{$name} = $self->{_das_segments}->{$key};
    } else { # Get fresh data
        my ($featref, $styleref, $segref) = $dasfact->fetch_all_Features( $slice, $type );
        $self->{_das_features}->{$key} = $featref;
        $self->{_das_styles}->{$key} = $styleref;
        $self->{_das_segments}->{$key} = $segref;
        $das_features{$name} = $featref;
        $das_styles{$name} = $styleref;
        $das_segments{$name} = $segref;
    }
  }

  return (\%das_features, \%das_styles, \%das_segments);
}

sub get_all_DASFeatures{
   my ($self, $source_type) = @_;


  if(!$self->adaptor()) {
    warning("Cannot retrieve features without attached adaptor");
    return [];
  }

  my %genomic_features = map { ( $_->adaptor->dsn => [ $_->fetch_all_Features($self, $source_type) ]  ) } $self->adaptor()->db()->_each_DASFeatureFactory;
  return \%genomic_features;

}

sub old_get_all_DASFeatures{
   my ($self,@args) = @_;

  if(!$self->adaptor()) {
    warning("Cannot retrieve features without attached adaptor");
    return [];
  }

  my %genomic_features =
      map { ( $_->adaptor->dsn => [ $_->fetch_all_by_Slice($self) ]  ) }
         $self->adaptor()->db()->_each_DASFeatureFactory;
  return \%genomic_features;

}


=head2 get_all_ExternalFeatures

  Arg [1]    : (optional) string $track_name
               If specified only features from ExternalFeatureAdaptors with
               the track name $track_name are retrieved.
               If not set, all features from every ExternalFeatureAdaptor are
               retrieved.
  Example    : @x_features = @{$slice->get_all_ExternalFeatures}
  Description: Retrieves features on this slice from external feature adaptors
  Returntype : listref of Bio::SeqFeatureI implementing objects in slice
               coordinates
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_ExternalFeatures {
   my ($self, $track_name) = @_;

   if(!$self->adaptor()) {
     warning("Cannot retrieve features without attached adaptor");
     return [];
   }

   my $features = [];

   my $xfa_hash = $self->adaptor->db->get_ExternalFeatureAdaptors;
   my @xf_adaptors = ();

   if($track_name) {
     #use a specific adaptor
     if(exists $xfa_hash->{$track_name}) {
       push @xf_adaptors, $xfa_hash->{$track_name};
     }
   } else {
     #use all of the adaptors
     push @xf_adaptors, values %$xfa_hash;
   }


   foreach my $xfa (@xf_adaptors) {
     push @$features, @{$xfa->fetch_all_by_Slice($self)};
   }

   return $features;
}


=head2 get_all_DitagFeatures

  Arg [1]    : (optional) string ditag type
  Arg [1]    : (optional) string logic_name
  Example    : @dna_dna_align_feats = @{$slice->get_all_DitagFeatures};
  Description: Retrieves the DitagFeatures of a specific type which overlap
               this slice with. If type is not defined, all features are
               retrieved.
  Returntype : listref of Bio::EnsEMBL::DitagFeatures
  Exceptions : warning if slice does not have attached adaptor
  Caller     : general
  Status     : Stable

=cut

sub get_all_DitagFeatures {
   my ($self, $type, $logic_name) = @_;

   if(!$self->adaptor()) {
     warning('Cannot get DitagFeatures without attached adaptor');
     return [];
   }

   my $dfa = $self->adaptor->db->get_DitagFeatureAdaptor();

   return $dfa->fetch_all_by_Slice($self, $type, $logic_name);
}




# GENERIC FEATURES (See DBAdaptor.pm)

=head2 get_generic_features

  Arg [1]    : (optional) List of names of generic feature types to return.
               If no feature names are given, all generic features are
               returned.
  Example    : my %features = %{$slice->get_generic_features()};
  Description: Gets generic features via the generic feature adaptors that
               have been added via DBAdaptor->add_GenricFeatureAdaptor (if
               any)
  Returntype : Hash of named features.
  Exceptions : none
  Caller     : none
  Status     : Stable

=cut

sub get_generic_features {

  my ($self, @names) = @_;

  if(!$self->adaptor()) {
    warning('Cannot retrieve features without attached adaptor');
    return [];
  }

  my $db = $self->adaptor()->db();

  my %features = ();   # this will hold the results

  # get the adaptors for each feature
  my %adaptors = %{$db->get_GenericFeatureAdaptors(@names)};

  foreach my $adaptor_name (keys(%adaptors)) {

    my $adaptor_obj = $adaptors{$adaptor_name};
    # get the features and add them to the hash
    my $features_ref = $adaptor_obj->fetch_all_by_Slice($self);
    # add each feature to the hash to be returned
    foreach my $feature (@$features_ref) {
      $features{$adaptor_name} = $feature;
    }
  }

  return \%features;

}
*/

/*
=head2 project_to_slice

  Arg [1]    : Slice to project to.
  Example    : my $chr_projection = $clone_slice->project_to_slice($chrom_slice);
                foreach my $segment ( @$chr_projection ){
                  $chr_slice = $segment->to_Slice();
                  print $clone_slice->seq_region_name(). ':'. $segment->from_start(). '-'.
                        $segment->from_end(). ' -> '.$chr_slice->seq_region_name(). ':'. $chr_slice->start().
                        '-'.$chr_slice->end().
                         $chr_slice->strand(). " length: ".($chr_slice->end()-$chr_slice->start()+1). "\n";
                }
  Description: Projection of slice to another specific slice. Needed for where we have multiple mappings
               and we want to state which one to project to.
  Returntype : list reference of Bio::EnsEMBL::ProjectionSegment objects which
               can also be used as [$start,$end,$slice] triplets.
  Exceptions : none
  Caller     : none
  Status     : At Risk

=cut
*/

Vector *Slice_projectToSlice(Slice *slice, Slice *toSlice) {

  if (!toSlice) {
    fprintf(stderr,"Slice argument is required\n");
    exit(1);
  }

  Vector *projection = Vector_new();
  SliceAdaptor *sa = (SliceAdaptor *)Slice_getAdaptor(slice);

  if (sa == NULL) {
    fprintf(stderr, "Warning: Cannot project without attached adaptor.\n");
    return projection;
  }

  if (Slice_getCoordSystem(slice) == NULL) {
    fprintf(stderr, "Warning: Cannot project without attached coord system.\n");
    return projection;
  }

  AssemblyMapperAdaptor *asma = DBAdaptor_getAssemblyMapperAdaptor(sa->dba);

  CoordSystem *cs = Slice_getCoordSystem(toSlice);
  CoordSystem *sliceCs = Slice_getCoordSystem(slice);
  IDType toSliceId = Slice_getSeqRegionId(toSlice);

  if (cs == NULL) {
    fprintf(stderr, "Warning: Cannot project without attached coord system.\n");
    return projection;
  }

  long currentStart = 1;

  // decompose this slice into its symlinked components.
  // this allows us to handle haplotypes and PARs
  Vector *normalSliceProj = SliceAdaptor_fetchNormalizedSliceProjection(sa, slice, 0);

  int i;
  for (i=0; i<Vector_getNumElement(normalSliceProj); i++) {
    ProjectionSegment *segment = Vector_getElementAt(normalSliceProj, i);

    Slice *normalSlice = ProjectionSegment_getToSlice(segment);

    sliceCs = Slice_getCoordSystem(normalSlice);

    AssemblyMapper *asmMapper = AssemblyMapperAdaptor_fetchByCoordSystems(asma, sliceCs, cs);

    // perform the mapping between this slice and the requested system
    MapperRangeSet *coords;

    if (asmMapper != NULL ) {
      coords = AssemblyMapper_map(asmMapper, Slice_getSeqRegionName(normalSlice),
                                             Slice_getStart(normalSlice),
                                             Slice_getEnd(normalSlice),
                                             Slice_getStrand(normalSlice),
                                             sliceCs, 0, toSlice);
    } else {
      coords = MapperRangeSet_new();
      MapperRangeSet_addRange(coords, (MapperRange *)MapperGap_new(Slice_getStart(normalSlice), Slice_getEnd(normalSlice), 0));
    }

    int lastRank = 0;
    // construct a projection from the mapping results and return it
    int j;
    for (j=0; j<MapperRangeSet_getNumRange(coords); j++) {
      MapperRange *coord = MapperRangeSet_getRangeAt(coords, j);

      long coordStart  = coord->start;
      long coordEnd    = coord->end;
      long length      = coordEnd - coordStart + 1;

      if (lastRank != ((MapperGap *)coord)->rank){
        currentStart = 1;
      }
      lastRank = ((MapperGap *)coord)->rank;

      // skip gaps
      if (coord->rangeType == MAPPERRANGE_COORD) {
        MapperCoordinate *mc = (MapperCoordinate *)coord;

        if (mc->id != toSliceId) { // for multiple mappings only get the correct one
          currentStart += length;
          continue;
        }

        CoordSystem *coordCs = mc->coordSystem;

        // If the normalised projection just ended up mapping to the
        // same coordinate system we were already in then we should just
        // return the original region.  This can happen for example, if we
        // were on a PAR region on Y which refered to X and a projection to
        // 'toplevel' was requested.
//        if($coord_cs->equals($slice_cs)) {
//          # trim off regions which are not defined
//          return $self->_constrain_to_region();
//        }

        // create slices for the mapped-to coord system
        Slice *mcSlice = SliceAdaptor_fetchBySeqRegionId(sa, mc->id, coordStart, coordEnd, mc->strand);

        long currentEnd = currentStart + length - 1;

        Vector_addElement(projection, ProjectionSegment_new(currentStart, currentEnd, mcSlice));
      }

      currentStart += length;
    }
  }


  // delete the cache as we may want to map to different set next time and old
  // results will be cached.
  //AssemblyMapperAdaptor_deleteCache(asma);

  return projection;
}


/*
=head2 get_all_synonyms

  Args [1]   : String external_db_name The name of the database to retrieve 
               the synonym for
  Args [2]   : Integer external_db_version The version of the database to retrieve 
               the synonym for. If not specified then we will ignore any versions
  Example    : my @alternative_names = @{$slice->get_all_synonyms()};
               @alternative_names = @{$slice->get_all_synonyms('EMBL')};
  Description: get a list of alternative names for this slice
  Returntype : reference to list of SeqRegionSynonym objects.
  Exception  : none
  Caller     : general
  Status     : At Risk

=cut

sub get_all_synonyms{
  my ($self, $external_db_name, $external_db_version) = @_;

  if ( !defined( $self->{'synonym'} ) ) {
    my $adap = $self->adaptor->db->get_SeqRegionSynonymAdaptor();
    $self->{'synonym'} =
      $adap->get_synonyms( $self->get_seq_region_id($self) );
  }
  
  if(! $external_db_name) {
    return $self->{'synonym'};
  }
  my @args =  ($external_db_version) ? 
              ($external_db_name, $external_db_version) : 
              ($external_db_name, undef, 1);
  my $external_db_id = $self->adaptor->db()->get_DBEntryAdaptor()->get_external_db_id(@args);
  if(!$external_db_id) {
    my $extra = ($external_db_version) ? "and version $external_db_version " : q{};
    throw "The external database $external_db_name ${extra}did not result in a valid identifier";
  }

  return [ grep { $_->external_db_id() == $external_db_id } @{$self->{synonym}} ];
}
*/

/*
=head2 add_synonym

  Args[0]    : synonym.
  Example    : $slice->add_synonym("alt_name");
  Description: add an alternative name for this slice
  Returntype : none
  Exception  : none
  Caller     : general
  Status     : At Risk

=cut

sub add_synonym{
  my $self = shift;
  my $syn = shift;
  my $external_db_id = shift;
  
  my $adap = $self->adaptor->db->get_SeqRegionSynonymAdaptor();
  if ( !defined( $self->{'synonym'} ) ) {
    $self->{'synonym'} = $self->get_all_synonyms();
  }
  my $new_syn = Bio::EnsEMBL::SeqRegionSynonym->new( #-adaptor => $adap,
                                                     -synonym => $syn,
                                                     -external_db_id => $external_db_id, 
                                                     -seq_region_id => $self->get_seq_region_id($self));

  push (@{$self->{'synonym'}}, $new_syn);

  return;
}
*/

/*
=head2 summary_as_hash

  Example       : $slice_summary = $slice->summary_as_hash();
  Description   : Retrieves a textual summary of this slice.
  Returns       : hashref of descriptive strings
=cut

sub summary_as_hash {
  my $self = shift;
  my %summary;
  $summary{'display_id'} = $self->display_id;
  $summary{'start'} = $self->start;
  $summary{'end'} = $self->end;
  $summary{'strand'} = $self->strand;
  $summary{'Is_circular'} = $self->is_circular ? "true" : "false";
  $summary{'region_name'} = $self->seq_region_name();
  return \%summary;
}

#
# Bioperl Bio::PrimarySeqI methods:
#

=head2 id

  Description: Included for Bio::PrimarySeqI interface compliance (0.7)

=cut

sub id { name(@_); }


=head2 display_id

  Description: Included for Bio::PrimarySeqI interface compliance (1.2)

=cut

sub display_id { name(@_); }


=head2 primary_id

  Description: Included for Bio::PrimarySeqI interface compliance (1.2)

=cut

sub primary_id { name(@_); }


=head2 desc

Description: Included for Bio::PrimarySeqI interface compliance (1.2)

=cut

sub desc{ return $_[0]->coord_system->name().' '.$_[0]->seq_region_name(); }


=head2 moltype

Description: Included for Bio::PrimarySeqI interface compliance (0.7)

=cut

sub moltype { return 'dna'; }

=head2 alphabet

  Description: Included for Bio::PrimarySeqI interface compliance (1.2)

=cut

sub alphabet { return 'dna'; }


=head2 accession_number

  Description: Included for Bio::PrimarySeqI interface compliance (1.2)

=cut

sub accession_number { name(@_); }


# sub DEPRECATED METHODS #
###############################################################################

=head1 DEPRECATED METHODS

=head2 get_all_AffyFeatures

  Description:  DEPRECATED, use functionality provided by the Ensembl
                Functional Genomics API instead.

=cut

sub get_all_AffyFeatures {
    deprecate( 'Use functionality provided by the '
        . 'Ensembl Functional Genomics API instead.' );
    throw('Can not delegate deprecated functionality.');

    # Old code:

#    my $self = shift;
#    my @arraynames = @_;
#
#    my $sa = $self->adaptor();
#    if ( ! $sa ) {
#        warning( "Cannot retrieve features without attached adaptor." );
#    }
#    my $fa = $sa->db()->get_AffyFeatureAdaptor();
#    my $features;
#
#    if ( @arraynames ) {
#        $features = $fa->fetch_all_by_Slice_arrayname( $self, @arraynames );
#    } else {
#        $features = $fa->fetch_all_by_Slice( $self );
#    }
#    return $features;
}

=head2 get_all_OligoFeatures

  Description:  DEPRECATED, use functionality provided by the Ensembl
                Functional Genomics API instead.

=cut

sub get_all_OligoFeatures {

    deprecate( 'Use functionality provided by the '
        . 'Ensembl Functional Genomics API instead.' );
    throw('Can not delegate deprecated functionality.');

    # Old code:

#    my $self = shift;
#    my @arraynames = @_;
#
#    my $sa = $self->adaptor();
#    if ( ! $sa ) {
#        warning( "Cannot retrieve features without attached adaptor." );
#    }
#    my $fa = $sa->db()->get_OligoFeatureAdaptor();
#    my $features;
#
#    if ( @arraynames ) {
#        $features = $fa->fetch_all_by_Slice_arrayname( $self, @arraynames );
#    } else {
#        $features = $fa->fetch_all_by_Slice( $self );
#    }
#    return $features;
}

=head2 get_all_OligoFeatures_by_type

  Description:  DEPRECATED, use functionality provided by the Ensembl
                Functional Genomics API instead.

=cut

sub get_all_OligoFeatures_by_type {

    deprecate( 'Use functionality provided by the '
        . 'Ensembl Functional Genomics API instead.' );
    throw('Can not delegate deprecated functionality.');

    # Old code:

#    my ($self, $type, $logic_name) = @_;
#
#    throw('Need type as parameter') if !$type;
#
#    my $sa = $self->adaptor();
#    if ( ! $sa ) {
#        warning( "Cannot retrieve features without attached adaptor." );
#    }
#    my $fa = $sa->db()->get_OligoFeatureAdaptor();
#
#    my $features = $fa->fetch_all_by_Slice_type( $self, $type, $logic_name );
#
#    return $features;
}

=head2 get_all_supercontig_Slices

  Description: DEPRECATED use get_tiling_path("NTcontig") instead

=cut


sub get_all_supercontig_Slices {
  my $self = shift;

  deprecate("Use get_tiling_path('NTcontig') instead");

  my $result = [];

  if( $self->adaptor() ) {
    my $superctg_names =
      $self->adaptor()->list_overlapping_supercontigs( $self );

    for my $name ( @$superctg_names ) {
      my $slice;
      $slice = $self->adaptor()->fetch_by_supercontig_name( $name );
      $slice->name( $name );
      push( @$result, $slice );
    }
  } else {
    warning( "Slice needs to be attached to a database to get supercontigs" );
  }

  return $result;
}





=head2 get_Chromosome

  Description: DEPRECATED use this instead:
               $slice_adp->fetch_by_region('chromosome',
                                           $slice->seq_region_name)

=cut

sub get_Chromosome {
  my $self = shift @_;

  deprecate("Use SliceAdaptor::fetch_by_region('chromosome'," .
            '$slice->seq_region_name) instead');

  my $csa = $self->adaptor->db->get_CoordSystemAdaptor();
  my ($top_cs) = @{$csa->fetch_all()};

  return $self->adaptor->fetch_by_region($top_cs->name(),
                                         $self->seq_region_name(),
                                         undef,undef,undef,
                                         $top_cs->version());
}



=head2 chr_name

  Description: DEPRECATED use seq_region_name() instead

=cut

sub chr_name{
  deprecate("Use seq_region_name() instead");
  seq_region_name(@_);
}



=head2 chr_start

  Description: DEPRECATED use start() instead

=cut

sub chr_start{
  deprecate('Use start() instead');
  start(@_);
}



=head2 chr_end

  Description: DEPRECATED use end() instead
  Returntype : int
  Exceptions : none
  Caller     : SliceAdaptor, general

=cut

sub chr_end{
  deprecate('Use end() instead');
  end(@_);
}
*/


/*
=head2 assembly_type

  Description: DEPRECATED use version instead

=cut
*/

char *Slice_getAssemblyType(Slice *slice) {
  fprintf(stderr, "Slice_getAssemblyType deprecated. Use CoordSystem_getVersion(Slice_getCoordSystem(slice)) instead.\n");
  return CoordSystem_getVersion(Slice_getCoordSystem(slice));
}


/*
=head2 get_tiling_path

  Description: DEPRECATED use project instead

=cut

sub get_tiling_path {
  my $self = shift;
  deprecate('Use $slice->project("seqlevel") instead.');
  return [];
}


=head2 dbID

  Description: DEPRECATED use SliceAdaptor::get_seq_region_id instead

=cut

sub dbID {
  my $self = shift;
  deprecate('Use SliceAdaptor::get_seq_region_id instead.');
  if(!$self->adaptor) {
    warning('Cannot retrieve seq_region_id without attached adaptor.');
    return 0;
  }
  return $self->adaptor->get_seq_region_id($self);
}


=head2 get_all_MapFrags

  Description: DEPRECATED use get_all_MiscFeatures instead

=cut

sub get_all_MapFrags {
  my $self = shift;
  deprecate('Use get_all_MiscFeatures instead');
  return $self->get_all_MiscFeatures(@_);
}

=head2 has_MapSet

  Description: DEPRECATED use get_all_MiscFeatures instead

=cut

sub has_MapSet {
  my( $self, $mapset_name ) = @_;
  deprecate('Use get_all_MiscFeatures instead');
  my $mfs = $self->get_all_MiscFeatures($mapset_name);
  return (@$mfs > 0);
}
*/

