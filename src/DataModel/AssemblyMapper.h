#ifndef __ASSEMBLYMAPPER_H__
#define __ASSEMBLYMAPPER_H__

#include "DataModelTypes.h"

#include "IDHash.h"
#include "Mapper.h"
#include "AdaptorTypes.h"
#include "MapperRangeSet.h"

#define MAXCHUNK 2000

struct AssemblyMapperStruct {
  char *type;
  AssemblyMapperAdaptor *adaptor;
  Mapper *mapper;
  IDHash *contigRegister;
  IDHash *chrChunkHash;
};

#define AssemblyMapper_setChrChunkHash(am,h) (am)->chrChunkHash = (h)
#define AssemblyMapper_getChrChunkHash(am) (am)->chrChunkHash

#define AssemblyMapper_setAdaptor(am, ad) (am)->adaptor = (ad)
#define AssemblyMapper_getAdaptor(am) (am)->adaptor

#define AssemblyMapper_setMapper(am, m) (am)->mapper = (m)
#define AssemblyMapper_getMapper(am) (am)->mapper

char *AssemblyMapper_setType(AssemblyMapper *am, char *type);
#define AssemblyMapper_getType(am) (am)->type

#define AssemblyMapper_setContigRegister(am, cr) (am)->contigRegister = (cr)
#define AssemblyMapper_getContigRegister(am) (am)->contigRegister


AssemblyMapper *AssemblyMapper_new(AssemblyMapperAdaptor *ama, char *type);
char *AssemblyMapper_setType(AssemblyMapper *am, char *type);
MapperRangeSet *AssemblyMapper_mapCoordinatesToAssembly(AssemblyMapper *am, IDType contigId,
                                                        int start, int end, int strand);
int AssemblyMapper_fastToAssembly(AssemblyMapper *am, IDType contigId,
                                  int start, int end, int strand, MapperCoordinate *retRange);
MapperRangeSet *AssemblyMapper_mapCoordinatesToRawContig(AssemblyMapper *am, IDType chrId,
                              int start, int end, int strand);
int AssemblyMapper_listContigIds(AssemblyMapper *am, IDType chrId, int start, int end, IDType **ids);
void AssemblyMapper_registerRegion(AssemblyMapper *am, IDType chrId, int start, int end);
int AssemblyMapper_registerRegionAroundContig(AssemblyMapper *am, IDType contigId, int left, int right);
int AssemblyMapper_haveRegisteredContig(AssemblyMapper *am, IDType id);
void AssemblyMapper_registerContig(AssemblyMapper *am, IDType id);
int AssemblyMapper_getChunkSize(AssemblyMapper *am);
void AssemblyMapper_chunkRegisterRegion(AssemblyMapper *am, IDType chrId,
                     int firstChunk, int lastChunk);

#endif
