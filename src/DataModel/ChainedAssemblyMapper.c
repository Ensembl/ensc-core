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

#define __CHAINEDASSEMBLYMAPPER_MAIN__
#include "ChainedAssemblyMapper.h"
#undef  __CHAINEDASSEMBLYMAPPER_MAIN__
#include "CoordPair.h"
#include "AssemblyMapperAdaptor.h"
/*
=head1 DESCRIPTION

The ChainedAssemblyMapper is an extension of the regular AssemblyMapper
that allows for mappings between coordinate systems that require
multi-step mapping.  For example if explicit mappings are defined
between the following coordinate systems,

  chromosome <-> contig
  contig     <-> clone

the ChainedAssemblyMapper would be able to perform implicit mapping
between the chromosome and clone coordinate systems.  This should be
transparent to the user of this module, and users should not even
realise that they are using a chained assembly mapper as opposed to a
normal assembly mapper.

=head1 METHODS

=cut
*/


/*
my $FIRST = 'first';
my $MIDDLE = 'middle';
my $LAST  = 'last';
*/

// 2^20 = approx 10^6
int CAM_CHUNKFACTOR = 20;

// max size of the pair cache in the mappers
int CAM_DEFAULT_MAX_PAIR_COUNT = 6000;

/*
=head2 new

  Arg [1]    : Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor
  Arg [2]    : Bio::EnsEMBL::CoordSystem $src_cs
  Arg [3]    : Bio::EnsEMBL::CoordSystem $int_cs
  Arg [4]    : Bio::EnsEMBL::CoordSystem $dst_cs
  Example    : Should use AssemblyMapperAdaptor->fetch_by_CoordSystems
  Description: Creates a new AssemblyMapper
  Returntype : Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor
  Exceptions : thrown if wrong number of coord_systems are provided
  Caller     : AssemblyMapperAdaptor
  Status     : Stable

=cut
*/

ChainedAssemblyMapper *ChainedAssemblyMapper_new(AssemblyMapperAdaptor *adaptor, Vector *coordSystems) {
  ChainedAssemblyMapper *cam;

  if ((cam = (ChainedAssemblyMapper *)calloc(1, sizeof(ChainedAssemblyMapper))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for ChainedAssemblyMapper\n");
    return NULL;
  }

  cam->objectType = CLASS_CHAINEDASSEMBLYMAPPER;

  cam->funcs = &chainedAssemblyMapperFuncs;

  Object_incRefCount(cam);

  ChainedAssemblyMapper_setAdaptor(cam, adaptor);

  AssemblyMapperAdaptor_cacheSeqIdsWithMultAssemblies(adaptor);

  // Set the component, intermediate and assembled coordinate systems
  if ( Vector_getNumElement(coordSystems) != 3 ) {
    fprintf(stderr, "Can only map between two coordinate systems %d were provided\n", Vector_getNumElement(coordSystems));
    exit(1);
  }

  ChainedAssemblyMapper_setFirstCoordSystem(cam, Vector_getElementAt(coordSystems, 0));
  ChainedAssemblyMapper_setMiddleCoordSystem(cam, Vector_getElementAt(coordSystems, 1));
  ChainedAssemblyMapper_setLastCoordSystem(cam, Vector_getElementAt(coordSystems, 2));

  // maps between first and intermediate coord systems
  Mapper *firstMidMapper = Mapper_new("first", 
                                      "middle",
                                      ChainedAssemblyMapper_getFirstCoordSystem(cam),
                                      ChainedAssemblyMapper_getMiddleCoordSystem(cam)
                                     );
  ChainedAssemblyMapper_setFirstMiddleMapper(cam, firstMidMapper);

  // maps between last and intermediate
  Mapper *lastMidMapper = Mapper_new("last", 
                                     "middle",
                                     ChainedAssemblyMapper_getLastCoordSystem(cam),
                                     ChainedAssemblyMapper_getMiddleCoordSystem(cam)
                                    );
  ChainedAssemblyMapper_setLastMiddleMapper(cam, lastMidMapper);

  // mapper that is actually used and is loaded by the mappings generated
  // by the other two mappers
  Mapper *firstLastMapper = Mapper_new("first", 
                                       "last", 
                                       ChainedAssemblyMapper_getFirstCoordSystem(cam),
                                       ChainedAssemblyMapper_getLastCoordSystem(cam)
                                      );
  ChainedAssemblyMapper_setFirstLastMapper(cam, firstLastMapper);

  // need registries to keep track of what regions are registered in source
  // and destination coordinate systems
  ChainedAssemblyMapper_setFirstRegistry(cam, RangeRegistry_new());
  ChainedAssemblyMapper_setLastRegistry(cam, RangeRegistry_new());

  ChainedAssemblyMapper_setMaxPairCount(cam, CAM_DEFAULT_MAX_PAIR_COUNT);

  return cam;
}




/*
=head2 register_all

  Arg [1]    : none
  Example    : $mapper->max_pair_count(10e6);
               $mapper->register_all();
  Description: Pre-registers all assembly information in this mapper.  The
               cache size should be set to a sufficiently large value
               so that all of the information can be stored.  This method
               is useful when *a lot* of mapping will be done in regions
               which are distributed around the genome.   After registration
               the mapper will consume a lot of memory but will not have to
               perform any SQL and will be faster.
  Returntype : none
  Exceptions : none
  Caller     : specialised programs doing a lot of mapping
  Status     : Stable

=cut
*/

void ChainedAssemblyMapper_registerAllImpl(ChainedAssemblyMapper *cam) {
  AssemblyMapperAdaptor *ama = ChainedAssemblyMapper_getAdaptor(cam);

  AssemblyMapperAdaptor_registerAllChained(ama, cam);

  return;
}



void ChainedAssemblyMapper_flushImpl(ChainedAssemblyMapper *cam) {

  RangeRegistry_flush( ChainedAssemblyMapper_getFirstRegistry(cam) );
  RangeRegistry_flush( ChainedAssemblyMapper_getLastRegistry(cam) );

  Mapper_flush( ChainedAssemblyMapper_getFirstMiddleMapper(cam) );
  Mapper_flush( ChainedAssemblyMapper_getLastMiddleMapper(cam) );
  Mapper_flush( ChainedAssemblyMapper_getFirstLastMapper(cam) );

  return;
}

/*
=head2 size

  Args       : none
  Example    : $num_of_pairs = $mapper->size();
  Description: return the number of pairs currently stored.
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/

int ChainedAssemblyMapper_getSize(ChainedAssemblyMapper *cam) {
  return ( Mapper_getPairCount( ChainedAssemblyMapper_getFirstLastMapper(cam) ) +
           Mapper_getPairCount( ChainedAssemblyMapper_getLastMiddleMapper(cam) ) +
           Mapper_getPairCount( ChainedAssemblyMapper_getFirstMiddleMapper(cam) ) );
}



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
  Arg [6]    : (optional) fastmap
  Arg [7]    : (optional) Bio::Ensembl::Slice
               The slice to transform TO
  Example    : @coords = $asm_mapper->map('X', 1_000_000, 2_000_000,
                                            1, $chr_cs);
  Description: Transforms coordinates from one coordinate system
               to another.
  Returntype : List of Bio::EnsEMBL::Mapper::Coordinate and/or
               Bio::EnsEMBL::Mapper:Gap objects
  Exceptions : thrown if the specified TO coordinat system is not one
               of the coordinate systems associated with this assembly mapper
  Caller     : general
  Status     : Stable

=cut
*/

/* Now in BaseAssemblyMapper
IDType ChainedAssemblyMapper_getSeqRegionId(ChainedAssemblyMapper *cam, char *seqRegionName, CoordSystem *cs) {
  // Not the most efficient thing to do making these temporary vectors to get one value, but hey its what the perl does!
  Vector *tmp = Vector_new();
  Vector_addElement(tmp, seqRegionName);

  AssemblyMapperAdaptor *adaptor = ChainedAssemblyMapper_getAdaptor(cam);

  Vector *idVec = AssemblyMapperAdaptor_seqRegionsToIds(adaptor, cs, tmp);

  IDType seqRegionId = *((IDType *)Vector_getElementAt(idVec, 0));

  Vector_free(tmp);
  Vector_free(idVec);
  // End of somewhat inefficient stuff

  return seqRegionId;
}
*/

MapperRangeSet *ChainedAssemblyMapper_mapImpl(ChainedAssemblyMapper *cam, char *frmSeqRegionName, long frmStart, long frmEnd, int frmStrand, 
                              CoordSystem *frmCs, int fastmap, Slice *toSlice) {
  Mapper *mapper       = ChainedAssemblyMapper_getFirstLastMapper(cam);
  CoordSystem *firstCs = ChainedAssemblyMapper_getFirstCoordSystem(cam);
  CoordSystem *lastCs  = ChainedAssemblyMapper_getLastCoordSystem(cam);

  int isInsert = (frmStart == frmEnd+1);

  char *frm;
  RangeRegistry *registry;

  IDType seqRegionId = ChainedAssemblyMapper_getSeqRegionId(cam, frmSeqRegionName, frmCs);

  // speed critical section:
  // try to do simple pointer equality comparisons of the coord system objects
  // first since this is likely to work most of the time and is much faster
  // than a function call

  if (frmCs == firstCs || (frmCs != lastCs && !CoordSystem_compare(frmCs, firstCs))) {
    frm = "first";
    registry = ChainedAssemblyMapper_getFirstRegistry(cam);
  } else if (frmCs == lastCs || !CoordSystem_compare(frmCs, lastCs)) {
    frm = "last";
    registry = ChainedAssemblyMapper_getLastRegistry(cam);
  } else {
    fprintf(stderr,"Coordinate system %s %s is neither the first nor the last coordinate system "
                   " of this ChainedAssemblyMapper\n", CoordSystem_getName(frmCs), CoordSystem_getVersion(frmCs) );
    exit(1);
  }

  // the minimum area we want to register if registration is necessary is
  // about 1MB. Break requested ranges into chunks of 1MB and then register
  // this larger region if we have a registry miss.

  // use bitwise shift for fast and easy integer multiplication and division
  long minStart, minEnd;

  if (isInsert) {
    minStart = ((frmEnd >> CAM_CHUNKFACTOR) << CAM_CHUNKFACTOR);
    minEnd   = (((frmStart >> CAM_CHUNKFACTOR) + 1) << CAM_CHUNKFACTOR) - 1 ;
  } else {
    minStart = ((frmStart >> CAM_CHUNKFACTOR) << CAM_CHUNKFACTOR);
    minEnd   = (((frmEnd >> CAM_CHUNKFACTOR) + 1) << CAM_CHUNKFACTOR) - 1 ;
  }

  // get a list of ranges in the requested region that have not been registered,
  // and register them at the same

  Vector *ranges;

  if (isInsert) {
    ranges = RangeRegistry_checkAndRegister(registry, seqRegionId, frmEnd, frmStart, minStart, minEnd, 1);
  } else {
    ranges = RangeRegistry_checkAndRegister(registry, seqRegionId, frmStart, frmEnd, minStart, minEnd, 1);
  }

  if (Vector_getNumElement(ranges)) {
    if (ChainedAssemblyMapper_getSize(cam) > ChainedAssemblyMapper_getMaxPairCount(cam)) {
      ChainedAssemblyMapper_flush(cam);

      if (isInsert) {
        ranges = RangeRegistry_checkAndRegister(registry, seqRegionId, frmEnd, frmStart, minStart, minEnd, 1);
      } else {
        ranges = RangeRegistry_checkAndRegister(registry, seqRegionId, frmStart, frmEnd, minStart, minEnd, 1);
      }
    }
    AssemblyMapperAdaptor *adaptor = ChainedAssemblyMapper_getAdaptor(cam);
    AssemblyMapperAdaptor_registerChained(adaptor, cam, frm, seqRegionId, ranges, toSlice);
  }

  MapperRangeSet *mrs;
  if (fastmap) {
    mrs = Mapper_fastMap(mapper, seqRegionId, frmStart, frmEnd, frmStrand, frm);
  } else {
    mrs = Mapper_mapCoordinates(mapper, seqRegionId, frmStart, frmEnd, frmStrand, frm);
  }

  // NIY: Need to tidy up elements too!!
  Vector_free(ranges);

  return mrs;
}


MapperRangeSet *ChainedAssemblyMapper_fastMapImpl(ChainedAssemblyMapper *cam, char *frmSeqRegionName, long frmStart, long frmEnd, int frmStrand, 
                              CoordSystem *frmCs, Slice *fakeSliceArg) {
  return ChainedAssemblyMapper_map(cam, frmSeqRegionName, frmStart, frmEnd, frmStrand, frmCs, 1, NULL);
}


/*
=head2 list_ids

  Arg [1]    : string $frm_seq_region
               The name of the sequence region of interest
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
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/


Vector *ChainedAssemblyMapper_listIdsImpl(ChainedAssemblyMapper *cam, char *frmSeqRegionName, long frmStart, long frmEnd, CoordSystem *frmCs) {

  int isInsert = (frmStart == frmEnd+1);

  //the minimum area we want to register if registration is necessary is
  //about 1MB. Break requested ranges into chunks of 1MB and then register
  //this larger region if we have a registry miss.

  //use bitwise shift for fast and easy integer multiplication and division
  long minStart, minEnd;

  if (isInsert) {
    minStart = ((frmEnd >> CAM_CHUNKFACTOR) << CAM_CHUNKFACTOR);
    minEnd   = (((frmStart >> CAM_CHUNKFACTOR) + 1) << CAM_CHUNKFACTOR) - 1;
  } else {
    minStart = ((frmStart >> CAM_CHUNKFACTOR) << CAM_CHUNKFACTOR);
    minEnd   = (((frmEnd >> CAM_CHUNKFACTOR) + 1) << CAM_CHUNKFACTOR) - 1;
  }

  IDType seqRegionId = ChainedAssemblyMapper_getSeqRegionId(cam, frmSeqRegionName, frmCs);


  if (!CoordSystem_compare(frmCs, ChainedAssemblyMapper_getFirstCoordSystem(cam))) {
    RangeRegistry *registry = ChainedAssemblyMapper_getFirstRegistry(cam);

    Vector *ranges;

    if (isInsert) {
      ranges = RangeRegistry_checkAndRegister(registry, seqRegionId, frmEnd, frmStart, minStart, minEnd, 1);
    } else {
      ranges = RangeRegistry_checkAndRegister(registry, seqRegionId, frmStart, frmEnd, minStart, minEnd, 1);
    }

    if (ranges) {
      AssemblyMapperAdaptor *adaptor = ChainedAssemblyMapper_getAdaptor(cam);
      AssemblyMapperAdaptor_registerChained(adaptor, cam, "first", seqRegionId, ranges, NULL);
    }

    MapperPairSet *mps = Mapper_listPairs( ChainedAssemblyMapper_getFirstLastMapper(cam), seqRegionId, frmStart, frmEnd, "first");

    return MapperPairSet_getToIds(mps) ;

  } else if (!CoordSystem_compare(frmCs, ChainedAssemblyMapper_getLastCoordSystem(cam))) {
    RangeRegistry *registry = ChainedAssemblyMapper_getLastRegistry(cam);

    Vector *ranges;

    if (isInsert) {
      ranges = RangeRegistry_checkAndRegister(registry, seqRegionId, frmEnd, frmStart, minStart, minEnd, 1);
    } else {
      ranges = RangeRegistry_checkAndRegister(registry, seqRegionId, frmStart, frmEnd, minStart, minEnd, 1);
    }

    if (Vector_getNumElement(ranges)) {
      AssemblyMapperAdaptor *adaptor = ChainedAssemblyMapper_getAdaptor(cam);
      AssemblyMapperAdaptor_registerChained(adaptor, cam, "last", seqRegionId, ranges, NULL);
    }

    MapperPairSet *mps = Mapper_listPairs( ChainedAssemblyMapper_getFirstLastMapper(cam), seqRegionId, frmStart, frmEnd, "last");

    return MapperPairSet_getFromIds(mps) ;
  } else {
    fprintf(stderr,"Coordinate system %s %s is neither the first nor the last coordinate system "
                    " of this ChainedAssemblyMapper\n", CoordSystem_getName(frmCs), CoordSystem_getVersion(frmCs) );
    exit(1);
  }
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
  Description: Retrieves a list of overlapping seq_region internal identifiers
               of another coordinate system.  This is the same as the
               list_ids method but uses seq_region names rather internal ids
  Returntype : List of strings
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut
*/

Vector *ChainedAssemblyMapper_listSeqRegionsImpl(ChainedAssemblyMapper *cam, char *frmSeqRegionName, long frmStart, long frmEnd, CoordSystem *frmCs) {

  //retrieve the seq_region names
  Vector *seqRegs = ChainedAssemblyMapper_listIds(cam, frmSeqRegionName, frmStart, frmEnd, frmCs);

  // The seq_regions are from the 'to' coordinate system not the
  // from coordinate system we used to obtain them
// SMJS toCs doesn't seem to be used
//  CoordSystem *toCs;
//  if (!CoordSystem_compare(frmCs, ChainedAssemblyMapper_getFirstCoordSystem(cam))) {
//    toCs = ChainedAssemblyMapper_getLastCoordSystem(cam);
//  } else {
//    toCs = ChainedAssemblyMapper_getFirstCoordSystem(cam);
//  }

  // convert them to names
  AssemblyMapperAdaptor *adaptor = ChainedAssemblyMapper_getAdaptor(cam);

  Vector *regions = AssemblyMapperAdaptor_seqIdsToRegions(adaptor, seqRegs);
  // Need to tidy up seqRegs;
  Vector_free(seqRegs);

  return regions;
}

void ChainedAssemblyMapper_freeImpl(ChainedAssemblyMapper *cam) {
  Object_errorUnimplementedMethod((Object*)cam, "ChainedAssemblyMapper_free");
}

/*
 Methods supplied to maintain polymorphism with AssemblyMapper there
 is no real assembled or component in the chained mapper, since the
 ordering is arbitrary and both ends might actually be assembled, but
 these methods provide convenient synonyms
*/

/*
=head2 mapper

  Args       : none
  Example    : $mapper = $cam->mapper();
  Description: return the first_last_mapper.
  Returntype : Bio::EnsEMBL::Mapper
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut
*/

/* Don't implement for now
sub mapper {
  my $self = shift;
  return $self->first_last_mapper();
}
*/

/*
=head2 map_coordinates_to_assembly

  Description: DEPRECATED, use map() instead.

=cut
*/
// Note changed to use contigName rather than contigId to match current _map functionality - need to change all calls to match this
MapperRangeSet *ChainedAssemblyMapper_mapCoordinatesToAssemblyImpl(ChainedAssemblyMapper *cam, char *contigName, long start, long end, int strand) {
  fprintf(stderr,"Deprecated method ChainedAssemblyMapper_mapCoordinatesToAssembly. Use ChainedAssemblyMapper_map instead\n");

  return ChainedAssemblyMapper_map(cam, contigName, start, end, strand, ChainedAssemblyMapper_getComponentCoordSystem(cam), 0, NULL);
}

/*
=head2 fast_to_assembly

  Description: DEPRECATED, use map() instead.

=cut
*/

// Note changed to use contigName rather than contigId to match current _map functionality - need to change all calls to match this
MapperRangeSet *ChainedAssemblyMapper_fastToAssemblyImpl(ChainedAssemblyMapper *cam, char *contigName, long start, long end, int strand) {

  fprintf(stderr,"Deprecated method ChainedAssemblyMapper_fastToAssembly. Use ChainedAssemblyMapper_map instead\n");

  return ChainedAssemblyMapper_map(cam, contigName, start, end, strand, ChainedAssemblyMapper_getComponentCoordSystem(cam), 0, NULL);
}

/*
=head2 map_coordinates_to_rawcontig

  Description: DEPRECATED, use map() instead.

=cut
*/

MapperRangeSet *ChainedAssemblyMapper_mapCoordinatesToRawContigImpl(ChainedAssemblyMapper *cam, char *chrName, long start, long end, int strand) {
  fprintf(stderr,"Deprecated method ChainedAssemblyMapper_mapCoordinatesToRawContig. Use ChainedAssemblyMapper_map instead\n");

// perl assumes first cs = assembled cs
  return ChainedAssemblyMapper_map(cam, chrName, start, end, strand, ChainedAssemblyMapper_getAssembledCoordSystem(cam), 0, NULL);
}

/*
=head2 list_contig_ids

  Description: DEPRECATED, use list_ids() instead.

=cut
*/
Vector *ChainedAssemblyMapper_listContigIdsImpl(ChainedAssemblyMapper *cam, char *chrName, long start, long end, int strand) {
  fprintf(stderr,"Deprecated method ChainedAssemblyMapper_listContigIds. Use ChainedAssemblyMapper_listIds instead\n");

  return ChainedAssemblyMapper_listIds(cam, chrName, start, end, ChainedAssemblyMapper_getAssembledCoordSystem(cam));
}

/*


=head2 map_coordinates_to_assembly

  DEPRECATED use map() instead

=cut

sub map_coordinates_to_assembly {
  my ($self, $contig_id, $start, $end, $strand) = @_;

  deprecate('Use map() instead.');

  #not sure if contig_id is seq_region_id or name...
  return $self->map($contig_id, $start, $end, $strand,
                   $self->contig_CoordSystem());

}


=head2 fast_to_assembly

  DEPRECATED use map() instead

=cut

sub fast_to_assembly {
  my ($self, $contig_id, $start, $end, $strand) = @_;

  deprecate('Use map() instead.');

  #not sure if contig_id is seq_region_id or name...
  return $self->map($contig_id, $start, $end, $strand,
                    $self->contig_CoordSystem());
}


=head2 map_coordinates_to_rawcontig

  DEPRECATED use map() instead

=cut

sub map_coordinates_to_rawcontig {
  my ($self, $chr_name, $start, $end, $strand) = @_;

  deprecate('Use map() instead.');

  return $self->map($chr_name, $start, $end, $strand,
                    $self->assembled_CoordSystem());

}

=head2 list_contig_ids
  DEPRECATED Use list_ids instead

=cut

sub list_contig_ids {
  my ($self, $chr_name, $start, $end) = @_;

  deprecate('Use list_ids() instead.');

  return $self->list_ids($chr_name, $start, $end,
                         $self->assembled_CoordSystem());
}
*/


/* Hopefully don't need these deprecated methods

=head2 in_assembly

  Deprecated. Use map() or list_ids() instead

=cut

sub in_assembly {
  my ($self, $object) = @_;

  deprecate('Use map() or list_ids() instead.');

  my $csa = $self->db->get_CoordSystemAdaptor();

  my $top_level = $csa->fetch_top_level();

  my $asma = $self->adaptor->fetch_by_CoordSystems($object->coord_system(),
                                                   $top_level);

  my @list = $asma->list_ids($object->seq_region(), $object->start(),
                             $object->end(), $object->coord_system());

  return (@list > 0);
}
*/
