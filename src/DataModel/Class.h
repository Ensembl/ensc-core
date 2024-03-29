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

#ifndef __CLASS_H__
#define __CLASS_H__

typedef struct ClassStruct Class;
typedef struct ClassHierarchyNodeStruct ClassHierarchyNode;
typedef enum ClassTypeEnum {
  CLASS_NONE,
  CLASS_OBJECT,
  CLASS_STATEMENTHANDLE,
  CLASS_MYSQLSTATEMENTHANDLE,
  CLASS_RESULTROW,
  CLASS_MYSQLRESULTROW,
  CLASS_SEQFEATURE,
  CLASS_EXON,
  CLASS_STICKYEXON,
  CLASS_TRANSCRIPT,
  CLASS_GENE,
  CLASS_FEATURESET,
  CLASS_SIMPLEFEATURE,
  CLASS_INTRON,
  CLASS_INTRONSUPPORTINGEVIDENCE,
  CLASS_REPEATFEATURE,
  CLASS_BASEALIGNFEATURE,
  CLASS_FEATUREPAIR,
  CLASS_DNADNAALIGNFEATURE,
  CLASS_DNAPEPALIGNFEATURE,
  CLASS_REPEATCONSENSUS,
  CLASS_ENSROOT,
  CLASS_DBENTRY,
  CLASS_ANALYSIS,
  CLASS_SLICE,
  CLASS_RAWCONTIG,
  CLASS_BASECONTIG,
  CLASS_SEQUENCE,
  CLASS_CHROMOSOME,
  CLASS_CLONE,
  CLASS_TRANSLATION,
  CLASS_PREDICTIONTRANSCRIPT,
  CLASS_ANNOTATEDSEQFEATURE,
  CLASS_DNAFRAG,
  CLASS_GENOMEDB,
  CLASS_GENOMICALIGN,
  CLASS_HOMOLOGY,
  CLASS_SYNTENYREGION,
  CLASS_VECTOR,
  CLASS_SPECIES,
  CLASS_BASEASSEMBLYMAPPER,
  CLASS_ASSEMBLYMAPPER,
  CLASS_CHAINEDASSEMBLYMAPPER,
  CLASS_TOPLEVELASSEMBLYMAPPER,
  CLASS_COORDSYSTEM,
  CLASS_ATTRIBUTE,
  CLASS_SEQEDIT,
  CLASS_PREDICTIONEXON,
  CLASS_NUMCLASS
} ClassType;

struct ClassHierarchyNodeStruct {
  int nSubClass;
  Class *mClass;
  ClassHierarchyNode **subClasses;
};

struct ClassStruct {
  ClassType type;
  char *name;
};


int Class_isDescendent(ClassType parentType, ClassType descType);
int Class_assertType(ClassType wantedType, ClassType actualType);
Class *Class_findByType(ClassType type);


#endif
