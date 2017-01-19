/*
 * Copyright [1999-2017] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

#ifndef __MAPPERRANGE_H__
#define __MAPPERRANGE_H__

typedef enum MapperRangeTypeEnum {
  MAPPERRANGE_NONE,
  MAPPERRANGE_GAP,
  MAPPERRANGE_COORD,
  MAPPERRANGE_INDEL
} MapperRangeType;

typedef struct MapperRangeStruct MapperRange;

#define MAPPERRANGE_DATA \
   int rangeType; \
   long start; \
   long end;

struct MapperRangeStruct {
  MAPPERRANGE_DATA
};


#define MapperRange_getLength(mr) ((mr)->end - (mr)->start + 1)

#endif
