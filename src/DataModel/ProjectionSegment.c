/*
 * Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
