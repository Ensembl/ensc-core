/*
 * Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
 * Copyright [2016-2017] EMBL-European Bioinformatics Institute
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
