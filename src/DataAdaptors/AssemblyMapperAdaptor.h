#ifndef __ASSEMBLYMAPPERADAPTOR_H__
#define __ASSEMBLYMAPPERADAPTOR_H__

#include "GenomicRange.h"
#include "AssemblyMapper.h"

#include "BaseAdaptor.h"
#include "AdaptorTypes.h"


struct AssemblyMapperAdaptorStruct {
  BASEADAPTOR_DATA
};

AssemblyMapperAdaptor *AssemblyMapperAdaptor_new(DBAdaptor *dba);

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
