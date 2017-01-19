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

#define __TOPLEVELASSEMBLYMAPPER_MAIN__
#include "TopLevelAssemblyMapper.h"
#undef  __TOPLEVELASSEMBLYMAPPER_MAIN__
#include "AssemblyMapperAdaptor.h"
#include "CoordSystemAdaptor.h"
#include "DBAdaptor.h"
#include "StrUtil.h"
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

  TopLevelAssemblyMapper *tlam = NULL;

  if (!topLevelCs || !otherCs) {
    return tlam;
  }

  if (!CoordSystem_getIsTopLevel(topLevelCs)) {
    fprintf(stderr,"%s is not the toplevel CoordSystem.\n", CoordSystem_getName(topLevelCs));
    exit(1);
  }
  if (CoordSystem_getIsTopLevel(otherCs)) {
    fprintf(stderr,"Other coordsystem (%s) should NOT be the toplevel CoordSystem.\n", CoordSystem_getName(otherCs));
    exit(1);
  }

  CoordSystemAdaptor *csAdaptor = DBAdaptor_getCoordSystemAdaptor(ama->dba);
  Vector *coordSystems          = CoordSystemAdaptor_fetchAll(csAdaptor);

  if ((tlam = (TopLevelAssemblyMapper *)calloc(1, sizeof(TopLevelAssemblyMapper))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for TopLevelAssemblyMapper\n");
    return NULL;
  }

  tlam->objectType = CLASS_TOPLEVELASSEMBLYMAPPER;

  tlam->funcs = &topLevelAssemblyMapperFuncs;

  Object_incRefCount(tlam);

  TopLevelAssemblyMapper_setAdaptor(tlam, ama);
  tlam->coordSystems = coordSystems;
  tlam->topLevelCs   = topLevelCs;
  tlam->otherCs      = otherCs;

  return tlam;
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


/* Now in BaseAssemblyMapper
IDType TopLevelAssemblyMapper_getSeqRegionId(TopLevelAssemblyMapper *tlam, char *seqRegionName, CoordSystem *cs) {
  // Not the most efficient thing to do making these temporary vectors to get one value, but hey its what the perl does!
  Vector *tmp = Vector_new();
  Vector_addElement(tmp, seqRegionName);

  AssemblyMapperAdaptor *adaptor = TopLevelAssemblyMapper_getAdaptor(tlam);

  Vector *idVec = AssemblyMapperAdaptor_seqRegionsToIds(adaptor, cs, tmp);

  IDType seqRegionId = *((IDType *)Vector_getElementAt(idVec, 0));

  Vector_free(tmp);
  Vector_free(idVec);
  // End of somewhat inefficient stuff

  return seqRegionId;
}
*/

MapperRangeSet *TopLevelAssemblyMapper_mapImpl(TopLevelAssemblyMapper *tlam, char *frmSeqRegionName, long frmStart, long frmEnd, int frmStrand, 
                                           CoordSystem *frmCs, int fastMap, Slice *fakeSliceArg) {

  if (CoordSystem_getIsTopLevel(frmCs)) {
    fprintf(stderr,"The toplevel CoordSystem can only be mapped TO, not FROM.\n");
  }

  IDType seqRegionId = TopLevelAssemblyMapper_getSeqRegionId(tlam, frmSeqRegionName, frmCs);

  // Not used in perl my $mapper      = $self->{'mapper'};
  // Not used in perl my $toplevel_cs = $self->{'toplevel_cs'};
  CoordSystem           *otherCs = TopLevelAssemblyMapper_getOtherCoordSystem(tlam);
  AssemblyMapperAdaptor *adaptor = TopLevelAssemblyMapper_getAdaptor(tlam);

  if (frmCs != otherCs && CoordSystem_compare(frmCs, otherCs)) {
    fprintf(stderr,"Coordinate system %s %s is neither the assembled nor the component coordinate system"
                   " of this AssemblyMapper\n", CoordSystem_getName(frmCs), CoordSystem_getVersion(frmCs));
  }

  Vector *coordSystems = TopLevelAssemblyMapper_getCoordSystems(tlam);
  
  CoordSystemAdaptor *csa = DBAdaptor_getCoordSystemAdaptor(adaptor->dba);

  //
  // TBD try to make this more efficient
  // 
  int fromRank = CoordSystem_getRank(otherCs);

  int i;
  for (i=0;i<Vector_getNumElement(coordSystems);i++) {
    CoordSystem *cs = Vector_getElementAt(coordSystems,i);

    if (CoordSystem_getRank(cs) >= fromRank) {
      break;
    }

    //check if a mapping path even exists to this coordinate system
    Vector *path = CoordSystemAdaptor_getMappingPath(csa, cs, otherCs);

    if (path && Vector_getNumElement(path) >0) {
      // Try to map to this coord system. If we get back any coordinates then
      // it is our 'toplevel' that we were looking for
      AssemblyMapper *mapper = AssemblyMapperAdaptor_fetchByCoordSystems(adaptor, otherCs, cs);

      if (fastMap) {
        MapperRangeSet *result = AssemblyMapper_fastMap(mapper, frmSeqRegionName, frmStart, frmEnd, frmStrand, frmCs, NULL);
  // NIY: Not sure what condition should be here?
        if (result && MapperRangeSet_getNumRange(result) > 0) {
          return result;
        }
        //Perl was: return @result if(@result);
      } else {
        MapperRangeSet *coords = AssemblyMapper_map(mapper, frmSeqRegionName, frmStart, frmEnd, frmStrand, frmCs, 0, NULL);

        if (MapperRangeSet_getNumRange(coords) > 1 || MapperRangeSet_getRangeAt(coords,0)->rangeType != MAPPERRANGE_GAP) {
          return coords;
        }
      }
    }
  }

  // the toplevel coordinate system for the region requested *is* the
  // requested region.
//  if ($fastmap) {
//    return ($seq_region_id,$frm_start, $frm_end, $frm_strand, $other_cs);
//  }
// NIY: Not sure what rank should be so set to 0
  MapperCoordinate *mc = MapperCoordinate_new(seqRegionId, frmStart, frmEnd, frmStrand, otherCs, 0 /*rank*/);

  MapperRangeSet *coords = MapperRangeSet_new();
  MapperRangeSet_addRange(coords, (MapperRange *)mc);

  return coords;
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

void TopLevelAssemblyMapper_flushImpl(TopLevelAssemblyMapper *tlam) {
// No op
}

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

MapperRangeSet *TopLevelAssemblyMapper_fastMapImpl(TopLevelAssemblyMapper *tlam, char *frmSeqRegionName, long frmStart, long frmEnd, int frmStrand, 
                                               CoordSystem *frmCs,Slice *fakeSliceArg) {
  return TopLevelAssemblyMapper_map(tlam, frmSeqRegionName, frmStart, frmEnd, frmStrand, frmCs, 1, fakeSliceArg);
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

CoordSystem *TopLevelAssemblyMapper_assembledCoordSystem(TopLevelAssemblyMapper *tlam) {
  return TopLevelAssemblyMapper_getTopLevelCoordSystem(tlam);
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

CoordSystem *TopLevelAssemblyMapper_componentCoordSystem(TopLevelAssemblyMapper *tlam) {
  return TopLevelAssemblyMapper_getOtherCoordSystem(tlam);
}


// private function which implements both list functions depending on whether seqRegions are passed in. As it returns a Vector this is possible in
// C as well, although its a bit odd
Vector * TopLevelAssemblyMapper_list(TopLevelAssemblyMapper *tlam, char *frmSeqRegionName, long frmStart, long frmEnd, 
                                     CoordSystem *frmCs, int seqRegionsFlag)  {

// Wasn't used in perl  my $mapper      = $self->{'mapper'};
//  Wasn't used in perl my $toplevel_cs = $self->{'toplevel_cs'};
  CoordSystem           *otherCs = TopLevelAssemblyMapper_getOtherCoordSystem(tlam);
  AssemblyMapperAdaptor *adaptor = TopLevelAssemblyMapper_getAdaptor(tlam);

  if (CoordSystem_getIsTopLevel(frmCs)) {
    fprintf(stderr, "The toplevel CoordSystem can only be mapped TO, not FROM.\n");
  }
  if (frmCs != otherCs && CoordSystem_compare(frmCs, otherCs)) {
    fprintf(stderr,"Coordinate system %s %s is neither the assembled nor the component coordinate system"
                   " of this AssemblyMapper\n", CoordSystem_getName(frmCs), CoordSystem_getVersion(frmCs));
  }

  Vector *coordSystems = TopLevelAssemblyMapper_getCoordSystems(tlam);
  
  CoordSystemAdaptor *csa = DBAdaptor_getCoordSystemAdaptor(adaptor->dba);

  //
  // TBD try to make this more efficient
  // 
  int fromRank = CoordSystem_getRank(otherCs);

  int i;
  for (i=0;i<Vector_getNumElement(coordSystems);i++) {
    CoordSystem *cs = Vector_getElementAt(coordSystems,i);

    if (CoordSystem_getRank(cs) >= fromRank) {
      break;
    }

    //check if a mapping path even exists to this coordinate system
    Vector *path = CoordSystemAdaptor_getMappingPath(csa, cs, otherCs);

// NIY: Not sure what condition should be here
    if (path && Vector_getNumElement(path) >0) {
      // Try to map to this coord system. If we get back any coordinates then
      // it is our 'toplevel' that we were looking for
      AssemblyMapper *mapper = AssemblyMapperAdaptor_fetchByCoordSystems(adaptor, otherCs, cs);

      Vector *result;

      IDType seqRegionId = TopLevelAssemblyMapper_getSeqRegionId(tlam, frmSeqRegionName, frmCs);

      if (seqRegionsFlag) {
        result = AssemblyMapper_listSeqRegions(mapper, frmSeqRegionName, frmStart, frmEnd, frmCs);
      } else {
        result = AssemblyMapper_listIds(mapper, frmSeqRegionName, frmStart, frmEnd, frmCs);
      }

// NIY: Not sure what condition should be
      if (result && Vector_getNumElement(result) > 0) {
        return result;
      }
      //Perl did: return @result if(@result);
    }
  }

  // the toplevel coordinate system for the region requested *is* the

// NIY: Odd name always returned no matter if ask for name or ids??????
  Vector *resVec = Vector_new();
  char *tmpStr;
  StrUtil_copyString(&tmpStr, frmSeqRegionName, 0);

  Vector_addElement(resVec, tmpStr);

  return resVec;

// NONE OF THE CODE BELOW HERE WILL EVER BE ACCESSED!!!
// NIY: No need to implement because return above means it was never accessed in perl

/*
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
*/
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

Vector *TopLevelAssemblyMapper_listSeqRegionsImpl(TopLevelAssemblyMapper *tlam, char *frmSeqRegionName, long frmStart, long frmEnd, CoordSystem *frmCs) {
  return TopLevelAssemblyMapper_list(tlam, frmSeqRegionName, frmStart, frmEnd, frmCs,1);
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

Vector *TopLevelAssemblyMapper_listIdsImpl(TopLevelAssemblyMapper *tlam, char *frmSeqRegionName, long frmStart, long frmEnd, CoordSystem *frmCs) {
  return TopLevelAssemblyMapper_list(tlam, frmSeqRegionName, frmStart, frmEnd, frmCs,0);
}


void TopLevelAssemblyMapper_freeImpl(TopLevelAssemblyMapper *tlam) {
  Object_errorUnimplementedMethod((Object*)tlam, "TopLevelAssemblyMapper_free");
}


void TopLevelAssemblyMapper_registerAllImpl(TopLevelAssemblyMapper *am) {
  Object_errorUnimplementedMethod((Object*)am, "TopLevelAssemblyMapper_registerAll");
}

MapperRangeSet *TopLevelAssemblyMapper_mapCoordinatesToAssemblyImpl(TopLevelAssemblyMapper *am, char *contigName, long start, long end, int strand) {
  Object_errorUnimplementedMethod((Object*)am, "TopLevelAssemblyMapper_mapCoordinatesToAssembly");
  return NULL;
}

MapperRangeSet *TopLevelAssemblyMapper_fastToAssemblyImpl(TopLevelAssemblyMapper *am, char *contigName, long start, long end, int strand) {
  Object_errorUnimplementedMethod((Object*)am, "TopLevelAssemblyMapper_fastToAssembly");
  return NULL;
}

MapperRangeSet *TopLevelAssemblyMapper_mapCoordinatesToRawContigImpl(TopLevelAssemblyMapper *am, char *chrName, long start, long end, int strand) {
  Object_errorUnimplementedMethod((Object*)am, "TopLevelAssemblyMapper_mapCoordinatesToRawContig");
  return NULL;
}

Vector *TopLevelAssemblyMapper_listContigIdsImpl(TopLevelAssemblyMapper *am, char *chrName, long start, long end, int strand) {
  Object_errorUnimplementedMethod((Object*)am, "TopLevelAssemblyMapper_listContigIds");
  return NULL;
}
