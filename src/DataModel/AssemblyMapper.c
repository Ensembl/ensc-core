#include "AssemblyMapper.h"
/*
=head1 DESCRIPTION

The AssemblyMapper is a database aware mapper which faciliates
conversion of coordinates between any two coordinate systems with an
relationship explicitly defined in the assembly table.  In the future
it may be possible to perform multiple step (implicit) mapping between
coordinate systems.

It is implemented using the Bio::EnsEMBL::Mapper object, which is a
generic mapper object between disjoint coordinate systems.

=head1 METHODS

=cut
*/


/*
my $ASSEMBLED = 'assembled';
my $COMPONENT = 'component';
*/

int DEFAULT_MAX_PAIR_COUNT = 1000;


/*
=head2 new

  Arg [1]    : Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor
  Arg [2]    : Bio::EnsEMBL::CoordSystem $asm_cs
  Arg [3]    : Bio::EnsEMBL::CoordSystem $cmp_cs
  Example    : Should use AssemblyMapperAdaptor->fetch_by_CoordSystems()
  Description: Creates a new AssemblyMapper
  Returntype : Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor
  Exceptions : Throws if multiple coord_systems are provided
  Caller     : AssemblyMapperAdaptor
  Status     : Stable

=cut
*/

AssemblyMapper *AssemblyMapper_new(AssemblyMapperAdaptor *adaptor, Vector *coordSystems) {
  AssemblyMapper *am;

  if ((am = (AssemblyMapper *)calloc(1, sizeof(AssemblyMapper))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for AssemblyMapper\n");
    return NULL;
  }

  AssemblyMapper_setAdaptor(am, adaptor);

  AssemblyMapperAdaptor_cacheSeqIdsWithMultAssemblies(adaptor);

  if ( Vector_getNumElement(coordSystems) != 2 ) {
    fprintf(stderr, "Can only map between two coordinate systems %d were provided\n", Vector_getNumElement(coordSystems));
    exit(1);
  }

  // Set the component and assembled coordinate systems
  AssemblyMapper_setAssembledCoordSystem(am, Vector_getElementAt(coordSystems, 0));
  AssemblyMapper_setComponentCoordSystem(am, Vector_getElementAt(coordSystems, 1));

  AssemblyMapper_setAssembledRegister(am, IDHash_new(IDHASH_MEDIUM));
  AssemblyMapper_setComponentRegister(am, IDHash_new(IDHASH_MEDIUM));

  // We load the mapper calling the 'ASSEMBLED' the 'from' coord system
  // and the 'COMPONENT' the 'to' coord system.
  AssemblyMapper_setMapper(am, Mapper_new("assembled", "component", AssemblyMapper_getAssembledCoordSystem(am), 
                                                                    AssemblyMapper_getComponentCoordSystem(am)));

  AssemblyMapper_setMaxPairCount(am, DEFAULT_MAX_PAIR_COUNT);

  return am;
}

/*
=head2 register_all

  Arg [1]    : None
  Example    : $mapper->max_pair_count(10e6);
               $mapper->register_all();
  Description: Pre-registers all assembly information in this
               mapper.  The cache size should be set to a
               sufficiently large value so that all of the
               information can be stored.  This method is useful
               when *a lot* of mapping will be done in regions
               which are distributed around the genome.  After
               registration the mapper will consume a lot of memory
               but will not have to perform any SQL and will be
               faster.
  Return type: None
  Exceptions : None
  Caller     : Specialised programs doing a lot of mapping.
  Status     : Stable

=cut
*/

void AssemblyMapper_registerAll(AssemblyMapper *am) {
  AssemblyMapperAdaptor *adaptor = AssemblyMapper_getAdaptor(am);

  AssemblyMapperAdaptor_registerAll(adaptor, am);
}

/*
=head2 map

  Arg [1]    : string $frm_seq_region
               The name of the sequence region to transform FROM.
  Arg [2]    : int $frm_start
               The start of the region to transform FROM.
  Arg [3]    : int $frm_end
               The end of the region to transform FROM.
  Arg [4]    : int $strand
               The strand of the region to transform FROM.
  Arg [5]    : Bio::EnsEMBL::CoordSystem
               The coordinate system to transform FROM
  Example    : @coords =
                $asm_mapper->map( 'X', 1_000_000, 2_000_000, 1,
                                  $chr_cs );
  Description: Transforms coordinates from one coordinate system to
               another.
  Return type: List of Bio::EnsEMBL::Mapper::Coordinate and/or
               Bio::EnsEMBL::Mapper:Gap objects.
  Exceptions : Throws if if the specified TO coordinat system is not
               one of the coordinate systems associated with this
               assembly mapper.
  Caller     : General
  Status     : Stable

=cut
*/

MapperRangeSet *AssemblyMapper_map(AssemblyMapper *am, char *frmSeqRegionName, long frmStart, long frmEnd, int frmStrand, CoordSystem *frmCs, Slice *toSlice) {
  Mapper *mapper  = AssemblyMapper_getMapper(am);
  CoordSystem *asmCs  = AssemblyMapper_getAssembledCoordSystem(am);
  CoordSystem *cmpCs  = AssemblyMapper_getComponentCoordSystem(am);
  AssemblyMapperAdaptor *adaptor = AssemblyMapper_getAdaptor(am);

  char *frm;

  IDType seqRegionId = AssemblyMapper_getSeqRegionId(am, frmSeqRegionName, frmCs);

  // Speed critical section:
  // Try to do simple pointer equality comparisons of the coord system
  // objects first since this is likely to work most of the time and is
  // much faster than a function call.

  if ( frmCs == cmpCs
       || ( frmCs != asmCs && !CoordSystem_compare(frmCs,  cmpCs)) ) {
    if ( !IDHash_contains(AssemblyMapper_getComponentRegister(am), seqRegionId) ) {
      AssemblyMapperAdaptor_registerComponent( adaptor, am, seqRegionId);
    }
    frm = "component";

  } else if ( frmCs == asmCs || !CoordSystem_compare(frmCs, asmCs) ) {
    // This can be probably be sped up some by only calling registered
    // assembled if needed.
    AssemblyMapperAdaptor_registerAssembled( adaptor, am, seqRegionId, frmStart, frmEnd);

    frm = "assembled";
  } else {
    fprintf(stderr,"Coordinate system %s %s is neither the assembled nor the component coordinate system of this AssemblyMapper\n",
            CoordSystem_getName(frmCs), CoordSystem_getVersion(frmCs) );

  }

  return Mapper_mapCoordinates( mapper, seqRegionId, frmStart, frmEnd, frmStrand, frm );
}


/*
=head2 flush

  Args       : None
  Example    : None
  Description: Remove all cached items from this AssemblyMapper.
  Return type: None
  Exceptions : None
  Caller     : AssemblyMapperAdaptor
  Status     : Stable

=cut
*/

void AssemblyMapper_flush(AssemblyMapper *am) {

  Mapper_flush( AssemblyMapper_getMapper(am) );

  // IDHash_free may not be correct function type
  IDHash_free( AssemblyMapper_getComponentRegister(am), IDHash_free);
  IDHash_free( AssemblyMapper_getAssembledRegister(am), IDHash_free);
}

/*
=head2 size

  Args       : None
  Example    : $num_of_pairs = $mapper->size();
  Description: Returns the number of pairs currently stored.
  Return type: int
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut
*/

int AssemblyMapper_getSize(AssemblyMapper *am) {

  return Mapper_getPairCount( AssemblyMapper_getMapper(am));
}

/*
=head2 fastmap

  Arg [1]    : string $frm_seq_region
               The name of the sequence region to transform FROM.
  Arg [2]    : int $frm_start
               The start of the region to transform FROM.
  Arg [3]    : int $frm_end
               The end of the region to transform FROM.
  Arg [4]    : int $strand
               The strand of the region to transform FROM.
  Arg [5]    : Bio::EnsEMBL::CoordSystem
               The coordinate system to transform FROM.
  Example    : @coords =
                $asm_mapper->map( 'X', 1_000_000, 2_000_000, 1,
                                  $chr_cs );
  Description: Transforms coordinates from one coordinate system to
               another.
  Return type: List of Bio::EnsEMBL::Mapper::Coordinate and/or
               Bio::EnsEMBL::Mapper:Gap objects.
  Exceptions : Throws if the specified TO coordinat system is not
               one of the coordinate systems associated with this
               assembly mapper.
  Caller     : General
  Status     : Stable

=cut
*/

MapperRangeSet *AssemblyMapper_fastMap(AssemblyMapper *am, char *frmSeqRegionName, long frmStart, long frmEnd, int frmStrand, CoordSystem *frmCs, Slice *toSlice) {
  Mapper *mapper  = AssemblyMapper_getMapper(am);
  CoordSystem *asmCs  = AssemblyMapper_getAssembledCoordSystem(am);
  CoordSystem *cmpCs  = AssemblyMapper_getComponentCoordSystem(am);
  AssemblyMapperAdaptor *adaptor = AssemblyMapper_getAdaptor(am);
  char *frm;

  IDType seqRegionId = AssemblyMapper_getSeqRegionId(am, frmSeqRegionName, frmCs);


  // Speed critical section:
  // Try to do simple pointer equality comparisons of the coord system
  // objects first since this is likely to work most of the time and is
  // much faster than a function call.

  if ( frmCs == cmpCs
       || ( frmCs != asmCs && !CoordSystem_compare(frmCs,  cmpCs)) ) {
    if ( !IDHash_contains(AssemblyMapper_getComponentRegister(am), seqRegionId) ) {
      AssemblyMapperAdaptor_registerComponent( adaptor, am, seqRegionId);
    }
    frm = "component";

  } else if ( frmCs == asmCs || !CoordSystem_compare(frmCs, asmCs) ) {
    // This can be probably be sped up some by only calling registered
    // assembled if needed.
    AssemblyMapperAdaptor_registerAssembled( adaptor, am, seqRegionId, frmStart, frmEnd);

    frm = "assembled";
  } else {
    fprintf(stderr,"Coordinate system %s %s is neither the assembled nor the component coordinate system of this AssemblyMapper\n",
            CoordSystem_getName(frmCs), CoordSystem_getVersion(frmCs) );

  }

  return Mapper_fastMap( mapper, seqRegionId, frmStart, frmEnd, frmStrand, frm );
}

/*
=head2 list_ids

  Arg [1]    : string $frm_seq_region
               The name of the sequence region of interest.
  Arg [2]    : int $frm_start
               The start of the region of interest.
  Arg [3]    : int $frm_end
               The end of the region to transform of interest.
  Arg [5]    : Bio::EnsEMBL::CoordSystem $frm_cs
               The coordinate system to obtain overlapping IDs of.
  Example    : foreach my $id (
                        $asm_mapper->list_ids( 'X', 1, 1000, $ctg_cs ) )
                { ... }
  Description: Retrieves a list of overlapping seq_region names of
               another coordinate system.  This is the same as the
               list_ids method but uses seq_region names rather
               internal IDs.
  Return type: List of strings.
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut
*/

IDType AssemblyMapper_getSeqRegionId(AssemblyMapper *am, char *seqRegionName, CoordSystem *cs) {
  // Not the most efficient thing to do making these temporary vectors to get one value, but hey its what the perl does!
  Vector *tmp = Vector_new();
  Vector_addElement(tmp, seqRegionName);

  AssemblyMapperAdaptor *adaptor = ChainedAssemblyMapper_getAdaptor(am);
  
  Vector *idVec = AssemblyMapperAdaptor_seqRegionsToIds(cs, tmp);

  IDType seqRegionId = *((IDType *)Vector_getElementAt(idVec, 0));

  Vector_free(tmp);
  Vector_free(idVec);
  // End of somewhat inefficient stuff

  return seqRegionId;
}

Vector *AssemblyMapper_listIds(AssemblyMapper *am, char *frmSeqRegionName, long frmStart, long frmEnd, CoordSystem *frmCs) {
  IDType seqRegionId = AssemblyMapper_getSeqRegionId(am, frmSeqRegionName, frmCs);

  if ( !CoordSystem_compare(frmCs, AssemblyMapper_getComponentCoordSystem(am) ) ) {

    if ( !AssemblyMapper_haveRegisteredComponent(am, seqRegionId) ) {
      AssemblyMapper_registeredComponent(am, seqRegionId);
    }

    // Pull out the 'from' identifiers of the mapper pairs.  The we
    // loaded the assembled side as the 'from' side in the constructor.

    MapperPairSet *mps = Mapper_listPairs( AssemblyMapper_getMapper(am), seqRegionId, frmStart, frmEnd, "component");

    return MapperPairSet_getFromIds(mps);
  } else if ( !CoordSystem_compare(frmCs, AssemblyMapper_getAssembledCoordSystem(am) ) ) {

    AssemblyMapper_registeredComponent(am, seqRegionId, frmStart, frmEnd);

    // Pull out the 'to' identifiers of the mapper pairs we loaded the
    // component side as the 'to' coord system in the constructor.

    MapperPairSet *mps = Mapper_listPairs( AssemblyMapper_getMapper(am), seqRegionId, frmStart, frmEnd, "assembled");

    return MapperPairSet_getToIds(mps);
  } else {

    fprintf(stderr, "Coordinate system %s %s is neither the assembled nor the component coordinate system of this AssemblyMapper\n",
                    CoordSystem_getName(frmCs), CoordSystem_getVersion(frmCs) );

  }
}



/*
=head2 list_seq_regions

  Arg [1]    : string $frm_seq_region
               The name of the sequence region of interest.
  Arg [2]    : int $frm_start
               The start of the region of interest.
  Arg [3]    : int $frm_end
               The end of the region to transform of interest.
  Arg [5]    : Bio::EnsEMBL::CoordSystem $frm_cs
               The coordinate system to obtain overlapping IDs of.
  Example    : foreach my $id (
                                 $asm_mapper->list_seq_regions(
                                                   'X', 1, 1000, $chr_cs
                                 ) ) { ... }
  Description: Retrieves a list of overlapping seq_region internal
               identifiers of another coordinate system.  This is
               the same as the list_seq_regions method but uses
               internal identfiers rather than seq_region strings.
  Return type: List of ints.
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut
*/

Vector *AssemblyMapper_listSeqRegions(AssemblyMapper *am, char *frmSeqRegionName, long frmStart, long frmEnd, CoordSystem *frmCs) {
  //retrieve the seq_region names
  Vector *seqRegs = AssemblyMapper_listIds(am, frmSeqRegionName, frmStart, frmEnd, frmCs);

// SMJS toCs doesn't seem to be used
//  CoordSystem *toCs;
//  if (!CoordSystem_compare(frmCs, AssemblyMapper_getAssembledCoordSystem(am))) {
//    toCs = AssemblyMapper_getComponentCoordSystem(am);
//  } else {
//    toCs = AssemblyMapper_getAssembledCoordSystem(am);
//  }

  // convert them to names
  AssemblyMapperAdaptor *adaptor = AssemblyMapper_getAdaptor(am);

  Vector *regions = AssemblyMapperAdaptor_seqIdsToRegions(adaptor, seqRegs);

  // Need to tidy up seqRegs;
  Vector_free(seqRegs);

  return regions;
}

/*
=head2 have_registered_component

  Arg [1]    : string $cmp_seq_region
               The name of the sequence region to check for
               registration.
  Example    : if ( $asm_mapper->have_registered_component('AL240214.1') ) {}
  Description: Returns true if a given component region has
               been registered with this assembly mapper.  This
               should only be called by this class or the
               AssemblyMapperAdaptor.  In other words, do not use
               this method unless you really know what you are
               doing.
  Return type: Boolean (0 or 1)
  Exceptions : Throws on incorrect arguments.
  Caller     : Internal, AssemblyMapperAdaptor
  Status     : Stable

=cut
*/

int AssemblyMapper_haveRegisteredComponent(AssemblyMapper *am, IDType cmpSeqRegionId) {
  IDHash *componentRegister = AssemblyMapper_getComponentRegister(am);

  if ( !IDHash_contains(componentRegister, cmpSeqRegionId) ) {
    return 0;
  }
  return 1;
}

/*
=head2 have_registered_assembled

  Arg [1]    : string $asm_seq_region
               The name of the sequence region to check for
               registration.
  Arg [2]    : int $chunk_id
               The chunk number of the provided seq_region to check
               for registration.
  Example    : if ( $asm_mapper->have_registered_component( 'X', 9 ) ) { }
  Description: Returns true if a given assembled region chunk
               has been registered with this assembly mapper.
               This should only be called by this class or the
               AssemblyMapperAdaptor.  In other words, do not use
               this method unless you really know what you are
               doing.
  Return type: Boolean (0 or 1)
  Exceptions : Throws on incorrect arguments
  Caller     : Internal, AssemblyMapperAdaptor
  Status     : Stable

=cut
*/

int AssemblyMapper_haveRegisteredAssembled(AssemblyMapper *am, IDType asmSeqRegionId, int chunkId) {
  IDHash *assembledRegister = AssemblyMapper_getAssembledRegister(am);

  if ( !IDHash_contains(assembledRegister, asmSeqRegionId) ) {
    return 0;
  }

  IDHash *chunkHash = IDHash_getValue(assembledRegister, asmSeqRegionId);

  if (!IDHash_contains(chunkHash, (IDType)chunkId)) {
    return 0;
  }

  return 1;
}


/*
=head2 register_component

  Arg [1]    : integer $cmp_seq_region
               The dbID of the component sequence region to
               register.
  Example    : $asm_mapper->register_component('AL312341.1');
  Description: Flags a given component sequence region as registered
               in this assembly mapper.  This should only be called
               by this class or the AssemblyMapperAdaptor.
  Return type: None
  Exceptions : Throws on incorrect arguments
  Caller     : Internal, AssemblyMapperAdaptor
  Status     : Stable

=cut
*/

void AssemblyMapper_registerComponent(AssemblyMapper *am, IDType cmpSeqRegionId) {
  IDHash *componentRegister = AssemblyMapper_getComponentRegister(am);

  if ( !IDHash_contains(componentRegister, cmpSeqRegionId) ) {
    IDHash_add(componentRegister, cmpSeqRegionId, &trueVal);
  }
}

/*
=head2 register_assembled

  Arg [1]    : integer $asm_seq_region
               The dbID of the sequence region to register.
  Arg [2]    : int $chunk_id
               The chunk number of the provided seq_region to register.
  Example    : $asm_mapper->register_assembled( 'X', 4 );
  Description: Flags a given assembled region as registered in this
               assembly mapper.  This should only be called by this
               class or the AssemblyMapperAdaptor.  Do not call this
               method unless you really know what you are doing.
  Return type: None
  Exceptions : Throws on incorrect arguments
  Caller     : Internal, AssemblyMapperAdaptor
  Status     : Stable

=cut
*/
void AssemblyMapper_registerAssembled(AssemblyMapper *am, IDType asmSeqRegionId, int chunkId) {
  IDHash *assembledRegister = AssemblyMapper_getAssembledRegister(am);

  if ( !IDHash_contains(assembledRegister, asmSeqRegionId) ) {
    IDHash_add(assembledRegister, asmSeqRegionId, IDHash_new(IDHASH_MEDIUM));
  }

  IDHash *chunkHash = IDHash_getValue(assembledRegister, asmSeqRegionId);

  if (!IDHash_contains(chunkHash, (IDType)chunkId)) {
    IDHash_add(chunkHash, (IDType)chunkId, &trueVal);
  }
}


/* Not implementing deprecated methods

=head2 in_assembly

  Description: DEPRECATED, use map() or list_ids() instead.

=cut

sub in_assembly {
  my ( $self, $object ) = @_;

  deprecate('Use map() or list_ids() instead.');

  my $csa = $self->db->get_CoordSystemAdaptor();

  my $top_level = $csa->fetch_top_level();

  my $asma =
    $self->adaptor->fetch_by_CoordSystems( $object->coord_system(),
                                           $top_level );

  my @list = $asma->list_ids( $object->seq_region(),
                              $object->start(),
                              $object->end(),
                              $object->coord_system() );

  return ( @list > 0 );
}

=head2 map_coordinates_to_assembly

  Description: DEPRECATED, use map() instead.

=cut

sub map_coordinates_to_assembly {
  my ( $self, $contig_id, $start, $end, $strand ) = @_;

  deprecate('Use map() instead.');

  # Not sure if contig_id is seq_region_id or name...
  return
    $self->map( $contig_id, $start, $end, $strand,
                $self->contig_CoordSystem() );

}

=head2 fast_to_assembly

  Description: DEPRECATED, use map() instead.

=cut

sub fast_to_assembly {
  my ( $self, $contig_id, $start, $end, $strand ) = @_;

  deprecate('Use map() instead.');

  # Not sure if contig_id is seq_region_id or name...
  return
    $self->map( $contig_id, $start, $end, $strand,
                $self->contig_CoordSystem() );
}

=head2 map_coordinates_to_rawcontig

  Description: DEPRECATED, use map() instead.

=cut

sub map_coordinates_to_rawcontig {
  my ( $self, $chr_name, $start, $end, $strand ) = @_;

  deprecate('Use map() instead.');

  return
    $self->map( $chr_name, $start, $end, $strand,
                $self->assembled_CoordSystem() );
}

=head2 list_contig_ids

  Description: DEPRECATED, use list_ids() instead.

=cut

sub list_contig_ids {
  my ( $self, $chr_name, $start, $end ) = @_;

  deprecate('Use list_ids() instead.');

  return
    $self->list_ids( $chr_name, $start, $end,
                     $self->assembled_CoordSystem() );
}

*/

