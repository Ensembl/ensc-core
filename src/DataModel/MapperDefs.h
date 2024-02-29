/*
 * Copyright [1999-2024] EMBL-European Bioinformatics Institute
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

#ifndef __MAPPERDEFS_H__
#define __MAPPERDEFS_H__

#define MAPPER_FROM_IND 0
#define MAPPER_TO_IND   1

typedef enum CoordSystemEnum {
  ASSEMBLY_COORDS,
  RAWCONTIG_COORDS,
  CDNA_COORDS,
  GENOMIC_COORDS
} CoordSystemType;

#endif
