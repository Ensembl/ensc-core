#include "TopLevelAssemblyMapper.h"
#include "AssemblyMapperAdaptor.h"
#include "CoordSystemAdaptor.h"
#include "DBAdaptor.h"
/*
=head1 DESCRIPTION

The TopLevelAssemblyMapper performs mapping between a provided
coordinate system and the toplevel pseudo cooordinate system.  The
toplevel coordinate system is not a real coordinate system, but
represents the highest coordinate system that can be mapped to in a
given region.  It is only possible to perform unidirectional mapping
using this mapper, because it does not make sense to map from the
toplevel coordinate system to another coordinate system.
*/

/*

=head2 new

  Arg [1]    : Bio::EnsEMBL::DBAdaptor $dbadaptor the adaptor for
               the database this mapper is using.
  Arg [2]    : Toplevel CoordSystem
  Arg [3]    : Other CoordSystem
  Description: Creates a new TopLevelAssemblyMapper object
  Returntype : Bio::EnsEMBL::DBSQL::TopLevelAssemblyMapper
  Exceptions : throws if any of the 3 arguments are missing/ not 
             : of the correct type.
  Caller     : Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor
  Status     : Stable

=cut
*/

TopLevelAssemblyMapper *TopLevelAssemblyMapper_new(AssemblyMapperAdaptor *ama, CoordSystem *topLevelCs, CoordSystem *otherCs) {
  if (!CoordSystem_getIsTopLevel(topLevelCs)) {
    fprintf(stderr,"%s is not the toplevel CoordSystem.\n", CoordSystem_getName(topLevelCs));
    exit(1);
  }
  if (!CoordSystem_getIsTopLevel(otherCs)) {
    fprintf(stderr,"Other coordsystem (%s) should NOT be the toplevel CoordSystem.\n", CoordSystem_getName(otherCs));
    exit(1);
  }

  CoordSystemAdaptor *csAdaptor = DBAdaptor_getCoordSystemAdaptor(ama->dba);
  Vector *coordSystems          = CoordSystemAdaptor_fetchAll(csAdaptor);

  TopLevelAssemblyMapper *tlam;
  if ((tlam = (TopLevelAssemblyMapper *)calloc(1, sizeof(TopLevelAssemblyMapper))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for TopLevelAssemblyMapper\n");
    return NULL;
  }

  tlam->coordSystems = coordSystems;
  tlam->topLevelCs   = topLevelCs;
  tlam->otherCs      = otherCs;
  TopLevelAssemblyMapper_setAdaptor(tlam, ama);

  return tlam;
}


#ifdef NIY
/*
=head2 map

  Arg [1]    : string $frm_seq_region
               The name of the sequence region to transform FROM
  Arg [2]    : int $frm_start
               The start of the region to transform FROM
  Arg [3]    : int $frm_end
               The end of the region to transform FROM
  Arg [4]    : int $strand
               The strand of the region to transform FROM
  Arg [5]    : Bio::EnsEMBL::CoordSystem
               The coordinate system to transform FROM
  Arg [6]    : if set will do a fastmap
  Example    : @coords = $mapper->map('X', 1_000_000, 2_000_000,
                                            1, $chr_cs);
  Description: Transforms coordinates from one coordinate system
               to another.
  Returntype : List of Bio::EnsEMBL::Mapper::Coordinate and/or
               Bio::EnsEMBL::Mapper:Gap objects
  Exceptions : thrown if if the specified TO coordinate system is not one
               of the coordinate systems associated with this mapper
  Caller     : general
  Status     : Stable

=cut
*/


sub map {
  throw('Incorrect number of arguments.') if(@_ != 6 && @_ != 7);

  my($self, $frm_seq_region_name, $frm_start, $frm_end, $frm_strand, $frm_cs,
    $fastmap) = @_;

  if($frm_cs->is_top_level()) {
    throw("The toplevel CoordSystem can only be mapped TO, not FROM.");
  }

  my @tmp;
  push @tmp, $frm_seq_region_name;
  my $seq_region_id = @{$self->adaptor()->seq_regions_to_ids($frm_cs, \@tmp)}[0];

  my $mapper      = $self->{'mapper'};
  my $toplevel_cs = $self->{'toplevel_cs'};
  my $other_cs    = $self->{'other_cs'};
  my $adaptor     = $self->adaptor;

  if($frm_cs != $other_cs && !$frm_cs->equals($other_cs)) {
    throw("Coordinate system " . $frm_cs->name . " " . $frm_cs->version .
          " is neither the assembled nor the component coordinate system " .
          " of this AssemblyMapper");
  }

  my $coord_systems = $self->{'coord_systems'};

  my $csa = $self->adaptor()->db()->get_CoordSystemAdaptor();

  //
  // TBD try to make this more efficient
  // 
  my $from_rank = $other_cs->rank();
  foreach my $cs (@$coord_systems) {
    last if($cs->rank >= $from_rank);

    //check if a mapping path even exists to this coordinate system
    my @mapping_path = @{ $csa->get_mapping_path( $cs, $other_cs ) };

    if(@mapping_path) {

      // Try to map to this coord system. If we get back any coordinates then
      // it is our 'toplevel' that we were looking for
      my $mapper = $adaptor->fetch_by_CoordSystems($other_cs, $cs);

      if($fastmap) {
        my @result = $mapper->fastmap($frm_seq_region_name, $frm_start, $frm_end,
                                      $frm_strand, $frm_cs);
        return @result if(@result);
      } else {
        my @coords = $mapper->map($frm_seq_region_name, $frm_start, $frm_end,
                                  $frm_strand, $frm_cs);

        if(@coords > 1 || !$coords[0]->isa('Bio::EnsEMBL::Mapper::Gap')) {
          return @coords;
        }
      }
    }
  }

  // the toplevel coordinate system for the region requested *is* the
  // requested region.
  if($fastmap) {
    return ($seq_region_id,$frm_start, $frm_end, $frm_strand, $other_cs);
  }
  return Bio::EnsEMBL::Mapper::Coordinate->new
    ($seq_region_id,$frm_start,$frm_end, $frm_strand, $other_cs);
}

//
// for polymorphism with AssemblyMapper
//
/*
=head2 flush

  Args       : none
  Example    : none
  Description: polymorphism with AssemblyMapper, does nothing
  Returntype : none
  Exceptions : none
  Status     : Stable

=cut
*/

sub flush {}

/*
=head2 fastmap

  Arg [1]    : string $frm_seq_region
               The name of the sequence region to transform FROM
  Arg [2]    : int $frm_start
               The start of the region to transform FROM
  Arg [3]    : int $frm_end
               The end of the region to transform FROM
  Arg [4]    : int $strand
               The strand of the region to transform FROM
  Arg [5]    : Bio::EnsEMBL::CoordSystem
               The coordinate system to transform FROM
  Example    : @coords = $mapper->fastmap('X', 1_000_000, 2_000_000,
                                            1, $chr_cs);
  Description: Transforms coordinates from one coordinate system
               to another.
  Returntype : List of Bio::EnsEMBL::Mapper::Coordinate and/or
               Bio::EnsEMBL::Mapper:Gap objects
  Exceptions : thrown if if the specified TO coordinate system is not one
               of the coordinate systems associated with this mapper
  Caller     : general
  Status     : Stable

=cut
*/

sub fastmap {
  my $self = shift;
  return $self->map(@_,1);
}

/*
=head2 assembled_CoordSystem

  Arg [1]    : none
  Example    : $cs = $mapper->assembled_CoordSystem
  Description: Retrieves the assembled CoordSystem from this mapper
  Returntype : Bio::EnsEMBL::CoordSystem
  Exceptions : none
  Caller     : internal, AssemblyMapperAdaptor
  Status     : Stable

=cut
*/

sub assembled_CoordSystem {
  my $self = shift;
  return $self->{'toplevel_cs'};
}

/*
=head2 component_CoordSystem

  Arg [1]    : none
  Example    : $cs = $mapper->component_CoordSystem
  Description: Retrieves the component CoordSystem from this  mapper
  Returntype : Bio::EnsEMBL::CoordSystem
  Exceptions : none
  Caller     : internal, AssemblyMapperAdaptor
  Status     : Stable

=cut
*/

sub component_CoordSystem {
  my $self = shift;
  return $self->{'other_cs'};
}


sub _list {
  my($self, $frm_seq_region_name, $frm_start, $frm_end, $frm_cs, $seq_regions) = @_;

  my $mapper      = $self->{'mapper'};
  my $toplevel_cs = $self->{'toplevel_cs'};
  my $other_cs    = $self->{'other_cs'};
  my $adaptor     = $self->adaptor;

  if($frm_cs->is_top_level()) {
    throw("The toplevel CoordSystem can only be mapped TO, not FROM.");
  }
  if($frm_cs != $other_cs && !$frm_cs->equals($other_cs)) {
    throw("Coordinate system " . $frm_cs->name . " " . $frm_cs->version .
          " is neither the assembled nor the component coordinate system " .
          " of this AssemblyMapper");
  }

  my $coord_systems = $self->{'coord_systems'};
  my $csa = $self->adaptor()->db()->get_CoordSystemAdaptor();

  //
  // TBD try to make this more efficient
  //
  my $from_rank = $other_cs->rank();
  foreach my $cs (@$coord_systems) {
    last if($cs->rank >= $from_rank);

    //check if a mapping path even exists to this coordinate system
    my @mapping_path = @{ $csa->get_mapping_path( $cs, $other_cs ) };

    if(@mapping_path) {

      // Try to map to this coord system. If we get back any coordinates then
      // it is our 'toplevel' that we were looking for
      my $mapper = $adaptor->fetch_by_CoordSystems($other_cs, $cs);

      my @result;

      my @tmp;
      push @tmp, $frm_seq_region_name;
      my $seq_region_id = @{$self->adaptor()->seq_regions_to_ids($frm_cs, \@tmp)}[0];

      if($seq_regions) {
        @result = $mapper->list_seq_regions($frm_seq_region_name, $frm_start,
                                            $frm_end, $frm_cs);
      } else {
        @result = $mapper->list_ids($frm_seq_region_name, $frm_start,
                                    $frm_end, $frm_cs);
      }

      return @result if(@result);
    }
  }

  // the toplevel coordinate system for the region requested *is* the

  return ($frm_seq_region_name);

// NONE OF THE CODE BELOW HERE WILL EVER BE ACCESSED!!!

  // requested region.
  if($seq_regions) {
    return ($frm_seq_region_name);
  }

  //this seems a bit silly and inefficient, but it is probably never
  //called anyway.
  my $slice_adaptor = $adaptor->db()->get_SliceAdaptor();
  my $slice = $slice_adaptor->fetch_by_region($other_cs->name(),
                                              $frm_seq_region_name,
                                              undef,undef,undef,$other_cs);
  return ($slice_adaptor->get_seq_region_id($slice));
}



/*
=head2 list_seq_regions

  Arg [1]    : string $frm_seq_region
               The name of the sequence region of interest
  Arg [2]    : int $frm_start
               The start of the region of interest
  Arg [3]    : int $frm_end
               The end of the region to transform of interest
  Arg [5]    : Bio::EnsEMBL::CoordSystem $frm_cs
               The coordinate system to obtain overlapping ids of
  Example    : foreach $id ($asm_mapper->list_ids('X',1,1000,$ctg_cs)) {...}
  Description: Retrieves a list of overlapping seq_region names
               of another coordinate system.  This is the same as the
               list_ids method but uses seq_region names rather internal ids
  Returntype : List of strings
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/

sub list_seq_regions {
  throw('Incorrect number of arguments.') if(@_ != 5);
  return _list(@_,1);
}


/*
=head2 list_ids

  Arg [1]    : string $frm_seq_region
               The name of the sequence region of interest.
  Arg [2]    : int $frm_start
               The start of the region of interest
  Arg [3]    : int $frm_end
               The end of the region to transform of interest
  Arg [5]    : Bio::EnsEMBL::CoordSystem $frm_cs
               The coordinate system to obtain overlapping ids of
  Example    : foreach $id ($asm_mapper->list_ids('X',1,1000,$chr_cs)) {...}
  Description: Retrieves a list of overlapping seq_region internal identifiers
               of another coordinate system.  This is the same as the
               list_seq_regions method but uses internal identfiers rather 
               than seq_region strings
  Returntype : List of ints
  Exceptions : thrown if the from CoordSystem is the toplevel coord system
               thrown if the from CoordSystem is not the one used in the mapper
  Caller     : general
  Status     : Stable

=cut
*/

sub list_ids {
  throw('Incorrect number of arguments.') if(@_ != 5);
  return _list(@_,0);
}


#endif

