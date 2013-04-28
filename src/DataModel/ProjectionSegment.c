#include "ProjectionSegment.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

ProjectionSegment *ProjectionSegment_new(long fromStart, long fromEnd, Slice *toSlice) {
  ProjectionSegment *segment;

  if ((segment = (ProjectionSegment *)calloc(1,sizeof(ProjectionSegment))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for ProjectionSegment\n");
    return NULL;
  }

  ProjectionSegment_setFromStart(segment, fromStart);
  ProjectionSegment_setFromEnd(segment, fromEnd);
  ProjectionSegment_setToSlice(segment, toSlice);

  return segment;
}

void ProjectionSegment_free(ProjectionSegment *segment) {
  free(segment);
}
