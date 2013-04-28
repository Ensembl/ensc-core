#ifndef __PROJECTIONSEGMENT_H__
#define __PROJECTIONSEGMENT_H__

#include "DataModelTypes.h"
#include "EnsC.h"
#include "Slice.h"

struct ProjectionSegmentStruct {
  Slice * toSlice;
  long    fromStart;
  long    fromEnd;
};

ProjectionSegment *ProjectionSegment_new();

#define ProjectionSegment_setFromStart(seg,s) (seg)->fromStart = (s)
#define ProjectionSegment_getFromStart(seg) (seg)->fromStart

#define ProjectionSegment_setFromEnd(seg,e) (seg)->fromEnd = (e)
#define ProjectionSegment_getFromEnd(seg) (seg)->fromEnd

#define ProjectionSegment_setToSlice(seg,slice) (seg)->toSlice = (slice)
#define ProjectionSegment_getToSlice(seg) (seg)->toSlice


void ProjectionSegment_free(ProjectionSegment *segment);


#endif
