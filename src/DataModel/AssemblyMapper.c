#include "AssemblyMapper.h"
#include "GenomicRange.h"
#include "Mapper.h"
#include "AssemblyMapperAdaptor.h"

char *DEFINED = "defined";

AssemblyMapper *AssemblyMapper_new(AssemblyMapperAdaptor *ama, char *type) {
  AssemblyMapper *am;
  Mapper *mapper;

  
  if ((am = (AssemblyMapper *)calloc(1,sizeof(AssemblyMapper))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for assemblymapper\n");
    return NULL;
  }
  
  AssemblyMapper_setContigRegister(am, IDHash_new(IDHASH_MEDIUM));
  AssemblyMapper_setChrChunkHash(am, IDHash_new(IDHASH_SMALL));
  AssemblyMapper_setAdaptor(am, ama);
  AssemblyMapper_setType(am, type);

  mapper = Mapper_new(RAWCONTIG_COORDS, ASSEMBLY_COORDS);
  AssemblyMapper_setMapper(am, mapper);
  return am;
}

char *AssemblyMapper_setType(AssemblyMapper *am, char *type) {
  if ((am->type = (char *)malloc(strlen(type)+1)) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for assmapper type\n");
    exit(1);
  }

  strcpy(am->type,type);

  return am->type;
}

MapperRangeSet *AssemblyMapper_mapCoordinatesToAssembly(AssemblyMapper *am, long contigId, 
                                                        int start, int end, int strand) {
  if( !IDHash_contains(AssemblyMapper_getContigRegister(am), contigId)) {
    AssemblyMapper_registerRegionAroundContig(am,contigId, 0, 0 );
  }

  return Mapper_mapCoordinates(AssemblyMapper_getMapper(am),contigId, start, 
                               end, strand, RAWCONTIG_COORDS);
}


int AssemblyMapper_fastToAssembly(AssemblyMapper *am, long contigId, 
                                  int start, int end, int strand, MapperCoordinate *retRange) {
  if (!IDHash_contains(AssemblyMapper_getContigRegister(am), contigId)) {
    AssemblyMapper_registerRegionAroundContig(am,contigId, 0, 0 );
  }

  return Mapper_fastMap(AssemblyMapper_getMapper(am),contigId, start, 
                        end, strand, RAWCONTIG_COORDS, retRange);
}


MapperRangeSet *AssemblyMapper_mapCoordinatesToRawcontig(AssemblyMapper *am, long chrId, 
                              int start, int end, int strand) {
  AssemblyMapper_registerRegion(am,chrId, start, end);
  
  return Mapper_mapCoordinates(AssemblyMapper_getMapper(am),chrId, start, 
                               end, strand, ASSEMBLY_COORDS);
}


int AssemblyMapper_listContigIds(AssemblyMapper *am, long chrId, int start, int end, long  **ids) {
  MapperPairSet *pairs;
  int nPair;
  int i;
  
  AssemblyMapper_registerRegion(am, chrId, start, end);
  
  pairs = Mapper_listPairs(AssemblyMapper_getMapper(am), chrId, 
                           start, end, ASSEMBLY_COORDS);
  
  if ((*ids = (long *)calloc(MapperPairSet_getNumPair(pairs),sizeof(long))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating ids array\n");
    exit(1);
  }
  
  nPair = MapperPairSet_getNumPair(pairs);
  for (i=0;i<MapperPairSet_getNumPair(pairs);i++) {
    (*ids)[i] = MapperPair_getUnit(MapperPairSet_getPairAt(pairs,i),MAPPER_FROM_IND)->id;
  }

  MapperPairSet_free(pairs,FALSE);

  return nPair;
}


void AssemblyMapper_registerRegion(AssemblyMapper *am, long chrId, int start, int end) {
  int firstChunk = (int)(start / AssemblyMapper_getChunkSize(am));
  int lastChunk  = (int)(end / AssemblyMapper_getChunkSize(am));
  
  AssemblyMapper_chunkRegisterRegion(am, chrId, firstChunk, lastChunk );
}



int AssemblyMapper_registerRegionAroundContig(AssemblyMapper *am, long contigId, int left, int right) {
  GenomicRange *range; 

  if (AssemblyMapper_haveRegisteredContig(am, contigId) && 
      left == 0 && right==0 ) {
    if (Mapper_listPairs(AssemblyMapper_getMapper(am), contigId, -1, -1, RAWCONTIG_COORDS)) {
      return 1;
    } else {
      return 0;
    }
  }
   
  range = AssemblyMapperAdaptor_registerContig(AssemblyMapper_getAdaptor(am),am, AssemblyMapper_getType(am),contigId);

  if(range) {
     AssemblyMapper_registerRegion(am, GenomicRange_getChrId(range),
                                       GenomicRange_getChrStart(range)-left,
                                       GenomicRange_getChrEnd(range)+right);
    GenomicRange_free(range);
    return 1;
  } else {
    return 0;
  }
}



int AssemblyMapper_haveRegisteredContig(AssemblyMapper *am, long id) {

  if( IDHash_contains(AssemblyMapper_getContigRegister(am),id)) {
    return 1;
  } else {
    return 0;
  }
}


void AssemblyMapper_registerContig(AssemblyMapper *am, long id) {
  IDHash_add(AssemblyMapper_getContigRegister(am),id,DEFINED);
}

int AssemblyMapper_getChunkSize(AssemblyMapper *am) {
  return 1000000;
}

void AssemblyMapper_chunkRegisterRegion(AssemblyMapper *am,long chrId,
                     int firstChunk, int lastChunk) {
  int i;
  IDHash *chunkHash = AssemblyMapper_getChrChunkHash(am);
  char *chunkFlags;

  if (lastChunk > MAXCHUNK) {
    fprintf(stderr, "ERROR: MAXCHUNK exceeded in chunkRegisterRegion\n");
    exit(1);
  }

  if (!IDHash_contains(chunkHash, chrId)) {
    if ((chunkFlags = (char *)calloc(MAXCHUNK, sizeof(char))) == NULL) {
      fprintf(stderr, "ERROR: Failed allocating chunkFlags\n");
      exit(1);
    }
    IDHash_add(chunkHash,chrId,chunkFlags);
  }


  chunkFlags = IDHash_getValue(chunkHash,chrId);

  for(i = firstChunk; i <= lastChunk; i++ ) {
    
    if (chunkFlags[i]) {
      continue;
    } else {
      int start = i * AssemblyMapper_getChunkSize(am);
      int end = start + AssemblyMapper_getChunkSize(am) - 1;
      chunkFlags[i] = 1;
      AssemblyMapperAdaptor_registerRegion(AssemblyMapper_getAdaptor(am),am,AssemblyMapper_getType(am),chrId,start,end);
    }
  }
}


/* Doesn't seem to be used

=head2 in_assembly

  Arg  1     : Bio::EnsEMBL::Clone or
               Bio::EnsEMBL::RawContig $object_in_assembly
                
  Example    : none
  Description: tests if the given Clone or RawContig object is in the 
               assembly
  Returntype : int 0,1
  Exceptions : argument type is checked
  Caller     : general

=cut


sub in_assembly {
  my ($self, $object) = @_;

  my @contigs;

  unless(ref $object) {
    $self->throw("$object is not an object reference");
  }

  if($object->isa("Bio::EnsEMBL::Clone")) {
    #get contigs from this clone
    @contigs = @{$object->get_all_Contigs()}; 
  } elsif ($object->isa("Bio::EnsEMBL::RawContig")) {
    #we already have the contig we need
    @contigs = ($object);
  } else {
    #object is not a clone or a raw contig
    $self->throw("$object is not a RawContig or Clone object");
  }

  #verify at least one of these contigs is mapped to the assembly
  foreach my $contig (@contigs) {
    if($self->register_region_around_contig( $contig->dbID(),
					     0, 0)) {
      return 1;
    }
  }

  #none of the contigs was in the assembly (golden path)
  return 0;
}
*/
