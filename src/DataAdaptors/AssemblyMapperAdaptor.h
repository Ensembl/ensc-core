#ifndef __ASSEMBLYMAPPERADAPTOR_H__
#define __ASSEMBLYMAPPERADAPTOR_H__

#include "GenomicRange.h"
#include "AssemblyMapper.h"

#include "BaseAdaptor.h"
#include "AdaptorTypes.h"

#include "StringHash.h"


struct AssemblyMapperAdaptorStruct {
  BASEADAPTOR_DATA
  StringHash *typeCache;
};

AssemblyMapperAdaptor *AssemblyMapperAdaptor_new(DBAdaptor *dba);
AssemblyMapper *AssemblyMapperAdaptor_fetchByType(AssemblyMapperAdaptor *ama, char *type);

void AssemblyMapperAdaptor_registerRegion(AssemblyMapperAdaptor *ama,
                                          AssemblyMapper *assMapper,
                                          char *assemblyType,
                                          long chrId,
                                          int chrStart,
                                          int chrEnd);

GenomicRange *AssemblyMapperAdaptor_registerContig(AssemblyMapperAdaptor *ama,
                                          AssemblyMapper *assMapper,
                                          char *assemblyType,
                                          long contigId);

#endif
