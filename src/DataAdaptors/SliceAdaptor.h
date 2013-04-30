#ifndef __SLICEADAPTOR_H__
#define __SLICEADAPTOR_H__

#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "Slice.h"
#include "StringHash.h"

struct SliceAdaptorStruct {
  BASEADAPTOR_DATA
  StringHash *srNameCache;
  IDHash *srIdCache;
  IDHash *asmExcCache;

};

#define SA_INCLUDE_NON_REFERENCE 0x2
#define SA_INCLUDE_DUPLICATES    0x4

SliceAdaptor *SliceAdaptor_new(DBAdaptor *dba);

Slice *SliceAdaptor_fetchByLocation(SliceAdaptor *sa, char *location, char *coordSystemName, char *coordSystemVersion, int noWarnings, int noFuzz);
void SliceAdaptor_parseLocationToValues(SliceAdaptor *sa, char *location, int noWarnings, int noErrors,
                                        char **retSeqRegionName, long *retStart, long *retEnd, int *retStrand );

Vector *SliceAdaptor_fetchNormalizedSliceProjection(SliceAdaptor *sa, Slice *slice, int filterProjections);
IDType SliceAdaptor_getSeqRegionId(SliceAdaptor *sa, Slice *slice);

void SliceAdaptor_buildExceptionCache(SliceAdaptor *sa);
Vector *SliceAdaptor_filterSliceProjections(SliceAdaptor *sa, Slice *slice, Vector *projections);

Slice *SliceAdaptor_fetchByFeature(SliceAdaptor *sa, SeqFeature *feature, int size, int isPercent);
Slice *SliceAdaptor_fetchBySeqRegionId(SliceAdaptor *sa, IDType seqRegionId, long start, long end, int strand);
Slice *SliceAdaptor_fetchByRegion(SliceAdaptor *sa, char *coordSystemName, char *inputSeqRegionName, long start, long end, int strand, char *version, int noFuzz);

Slice *SliceAdaptor_fetchByName(SliceAdaptor *sa, char *name);

Vector *SliceAdaptor_fetchAll(SliceAdaptor *sa, char *csName, char *csVersion, int flags);
Slice *SliceAdaptor_fetchByTopLevelLocation(SliceAdaptor *sa, char *location, int noWarnings, int noFuzz);




#endif
