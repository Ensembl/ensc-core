#define __BASEASSEMBLYMAPPER_MAIN__
#include "BaseAssemblyMapper.h"
#undef __BASEASSEMBLYMAPPER_MAIN__
#include "AssemblyMapperAdaptor.h"

/*
Abstract base class for AssemblyMapper classes - very limited functionality (just utility method getSeqRegionId)
*/

IDType BaseAssemblyMapper_getSeqRegionId(BaseAssemblyMapper *am, char *seqRegionName, CoordSystem *cs) {
  // Not the most efficient thing to do making these temporary vectors to get one value, but hey its what the perl does!
  Vector *tmp = Vector_new();
  Vector_addElement(tmp, seqRegionName);

  AssemblyMapperAdaptor *adaptor = AssemblyMapper_getAdaptor(am);
  
  Vector *idVec = AssemblyMapperAdaptor_seqRegionsToIds(adaptor, cs, tmp);

  IDType seqRegionId = *((IDType *)Vector_getElementAt(idVec, 0));

  Vector_free(tmp);
  Vector_free(idVec);
  // End of somewhat inefficient stuff

  return seqRegionId;
}

void BaseAssemblyMapper_registerAllImpl(BaseAssemblyMapper *am) {
  Object_errorUnimplementedMethod(am, "BaseAssemblyMapper_registerAll");
}

void BaseAssemblyMapper_freeImpl(BaseAssemblyMapper *am) {
  Object_errorUnimplementedMethod(am, "BaseAssemblyMapper_free");
}

MapperRangeSet *BaseAssemblyMapper_mapImpl(BaseAssemblyMapper *am, char *frmSeqRegionName, long frmStart, long frmEnd, int frmStrand, 
                                           CoordSystem *frmCs, int fakeFastMapFlag, Slice *toSlice) {
  Object_errorUnimplementedMethod(am, "BaseAssemblyMapper_map");
  return NULL;
}

void BaseAssemblyMapper_flushImpl(BaseAssemblyMapper *am) {
  Object_errorUnimplementedMethod(am, "BaseAssemblyMapper_flush");
}

MapperRangeSet *BaseAssemblyMapper_fastMapImpl(BaseAssemblyMapper *am, char *frmSeqRegionName, long frmStart, long frmEnd, int frmStrand, 
                                               CoordSystem *frmCs, Slice *toSlice) {
  Object_errorUnimplementedMethod(am, "BaseAssemblyMapper_fastMap");
  return NULL;
}

Vector *BaseAssemblyMapper_listIdsImpl(BaseAssemblyMapper *am, char *frmSeqRegionName, long frmStart, long frmEnd, CoordSystem *frmCs) {
  Object_errorUnimplementedMethod(am, "BaseAssemblyMapper_listIds");
  return NULL;
}

Vector *BaseAssemblyMapper_listSeqRegionsImpl(BaseAssemblyMapper *am, char *frmSeqRegionName, long frmStart, long frmEnd, CoordSystem *frmCs) {
  Object_errorUnimplementedMethod(am, "BaseAssemblyMapper_listSeqRegions");
  return NULL;
}

MapperRangeSet *BaseAssemblyMapper_mapCoordinatesToAssemblyImpl(BaseAssemblyMapper *am, char *contigName, long start, long end, int strand) {
  Object_errorUnimplementedMethod(am, "BaseAssemblyMapper_mapCoordinatesToAssembly");
  return NULL;
}

MapperRangeSet *BaseAssemblyMapper_fastToAssemblyImpl(BaseAssemblyMapper *am, char *contigName, long start, long end, int strand) {
  Object_errorUnimplementedMethod(am, "BaseAssemblyMapper_fastToAssembly");
  return NULL;
}

MapperRangeSet *BaseAssemblyMapper_mapCoordinatesToRawContigImpl(BaseAssemblyMapper *am, char *chrName, long start, long end, int strand) {
  Object_errorUnimplementedMethod(am, "BaseAssemblyMapper_mapCoordinatesToRawContig");
  return NULL;
}

Vector *BaseAssemblyMapper_listContigIdsImpl(BaseAssemblyMapper *am, char *chrName, long start, long end, int strand) {
  Object_errorUnimplementedMethod(am, "BaseAssemblyMapper_listContigIds");
  return NULL;
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
*/

