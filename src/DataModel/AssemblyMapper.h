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

#define AssemblyMapper_setChrChunkHash(am,h) (am)->chrChunkHash = h
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
MapperRangeSet *AssemblyMapper_mapCoordinatesToAssembly(AssemblyMapper *am, int64 contigId,
                                                        int start, int end, int strand);
int AssemblyMapper_fastToAssembly(AssemblyMapper *am, int64 contigId,
                                  int start, int end, int strand, MapperCoordinate *retRange);
MapperRangeSet *AssemblyMapper_mapCoordinatesToRawContig(AssemblyMapper *am, int64 chrId,
                              int start, int end, int strand);
int AssemblyMapper_listContigIds(AssemblyMapper *am, int64 chrId, int start, int end, int64 **ids);
void AssemblyMapper_registerRegion(AssemblyMapper *am, int64 chrId, int start, int end);
int AssemblyMapper_registerRegionAroundContig(AssemblyMapper *am, int64 contigId, int left, int right);
int AssemblyMapper_haveRegisteredContig(AssemblyMapper *am, int64 id);
void AssemblyMapper_registerContig(AssemblyMapper *am, int64 id);
int AssemblyMapper_getChunkSize(AssemblyMapper *am);
void AssemblyMapper_chunkRegisterRegion(AssemblyMapper *am, int64 chrId,
                     int firstChunk, int lastChunk);











#endif
