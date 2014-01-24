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

#ifndef __BASEADAPTOR_H__
#define __BASEADAPTOR_H__

#include <mysql.h>
#include <string.h>
#include <stdlib.h>

#include "AdaptorTypes.h"
//#include "DBAdaptor.h"
#include "StatementHandle.h"
#include "DBAdaptor.h"
#include "Vector.h"
#include "Slice.h"
#include "AssemblyMapper.h"


#ifdef __MAIN_C__
char *Adaptor_TypeStrings[] = {
  "NONE",
  "GENE",
  "TRANSCRIPT",
  "TRANSLATION",
  "EXON",
  "ANALYSIS",
  "CLONE",
  "RAWCONTIG",
  "SLICE",
  "CHROMOSOME",
  "COORDSYSTEM",
  "METACONTAINER",
  "SEQUENCE",
  "DNAALIGNFEATURE",
  "ASSEMBLYMAPPER",
  "SIMPLEFEATURE",
  "PROTEINALIGNFEATURE",
  "REPEATFEATURE",
  "REPEATCONSENSUS",
  "PREDICTIONTRANSCRIPT",
  "SUPPORTINGFEATURE",
  "DBENTRY",
  "GENOMEDB",
  "GENOMICALIGN",
  "DNAFRAG",
  "COMPARADNAALIGNFEATURE",
  "SYNTENY",
  "SYNTENYREGION",
  "HOMOLOGY",
  "METACOORDCONTAINER",
  "TRANSCRIPTSUPPPORTINGFEATURE",
  "ATTRIBUTE",
  "INTRONSUPPORTINGEVIDENCE",
  "PREDICTIONEXON",
  "CACHINGSEQUENCE"
};

#else
 extern char *Adaptor_TypeStrings[];
#endif

enum Adaptor_Types {
  NONE,
  GENE_ADAPTOR, 
  TRANSCRIPT_ADAPTOR, 
  TRANSLATION_ADAPTOR, 
  EXON_ADAPTOR,
  ANALYSIS_ADAPTOR,
  CLONE_ADAPTOR,
  RAWCONTIG_ADAPTOR,
  SLICE_ADAPTOR,
  CHROMOSOME_ADAPTOR,
  COORDSYSTEM_ADAPTOR,
  META_CONTAINER,
  SEQUENCE_ADAPTOR,
  DNAALIGNFEATURE_ADAPTOR,
  ASSEMBLYMAPPER_ADAPTOR,
  SIMPLEFEATURE_ADAPTOR,
  PROTEINALIGNFEATURE_ADAPTOR,
  REPEATFEATURE_ADAPTOR,
  REPEATCONSENSUS_ADAPTOR,
  PREDICTIONTRANSCRIPT_ADAPTOR,
  SUPPORTINGFEATURE_ADAPTOR,
  DBENTRY_ADAPTOR,
  GENOMEDB_ADAPTOR,
  GENOMICALIGN_ADAPTOR,
  DNAFRAG_ADAPTOR,
  COMPARADNAALIGNFEATURE_ADAPTOR,
  SYNTENY_ADAPTOR,
  SYNTENYREGION_ADAPTOR,
  HOMOLOGY_ADAPTOR,
  METACOORD_CONTAINER,
  TRANSCRIPTSUPPORTINGFEATURE_ADAPTOR,
  ATTRIBUTE_ADAPTOR,
  INTRONSUPPORTINGEVIDENCE_ADAPTOR,
  PREDICTIONEXON_ADAPTOR,
  CACHINGSEQUENCE_ADAPTOR,
};

typedef char * NameTableType[][2];

typedef int      (*BaseAdaptor_StoreFunc)(BaseAdaptor *bfa, Vector *features);
typedef NameTableType *(*BaseAdaptor_GetTablesFunc)(void);
typedef char **  (*BaseAdaptor_GetColumnsFunc)(void);
typedef Vector *    (*BaseAdaptor_ObjectsFromStatementHandleFunc)(BaseAdaptor *bfa,StatementHandle *sth,
                                                                      AssemblyMapper *mapper, Slice *slice);
typedef char *   (*BaseAdaptor_FinalClauseFunc)(void);
typedef char *   (*BaseAdaptor_DefaultWhereClauseFunc)(void);
typedef NameTableType *(*BaseAdaptor_LeftJoinFunc)(void);

typedef StatementHandle *(*BaseAdaptor_PrepareFunc)(BaseAdaptor *ba, char *qStr,size_t len);

#define NAME 0
#define SYN  1

#define BASEADAPTOR_DEF(DBADAPTOR_TYPE) \
  DBADAPTOR_TYPE *dba; \
  int adaptorType; \
  IDType speciesId; \
  int straightJoinFlag; \
  BaseAdaptor_PrepareFunc prepare; \
  BaseAdaptor_StoreFunc store; \
  BaseAdaptor_GetTablesFunc getTables; \
  BaseAdaptor_GetColumnsFunc getColumns; \
  BaseAdaptor_ObjectsFromStatementHandleFunc objectsFromStatementHandle; \
  BaseAdaptor_FinalClauseFunc finalClause; \
  BaseAdaptor_LeftJoinFunc leftJoin; \
  BaseAdaptor_DefaultWhereClauseFunc defaultWhereClause;

#define BASEADAPTOR_DATA BASEADAPTOR_DEF(DBAdaptor)

struct BaseAdaptorStruct {
  BASEADAPTOR_DATA
};

void BaseAdaptor_init(BaseAdaptor *ba, DBAdaptor *dba, int adaptorType);
StatementHandle *BaseAdaptor_prepare(BaseAdaptor *ba, char *qStr,size_t len);
NameTableType *BaseAdaptor_getTables(void);
char **BaseAdaptor_getColumns(void);
char *BaseAdaptor_finalClause(void);
NameTableType *BaseAdaptor_leftJoin(void);
char *BaseAdaptor_defaultWhereClause(void);
int BaseAdaptor_store(BaseAdaptor *bfa, Vector *features);
Vector *BaseAdaptor_objectsFromStatementHandle(BaseAdaptor *bfa, StatementHandle *sth,AssemblyMapper *mapper, Slice *slice);

Vector *BaseAdaptor_genericFetch(BaseAdaptor *ba, char *constraint, AssemblyMapper *mapper, Slice *slice);
Vector *BaseAdaptor_listDbIDs(BaseAdaptor *ba, char *table, char *pk, int ordered);
int BaseAdaptor_genericCount(BaseAdaptor *ba, char *constraint);
void BaseAdaptor_generateSql(BaseAdaptor *ba, char *constraint, char **inputColumns, char *sql);
SeqFeature *BaseAdaptor_fetchByDbID(BaseAdaptor *ba, IDType id);
Vector *BaseAdaptor_fetchAllByDbIDList(BaseAdaptor *ba, Vector *idList, Slice *slice);
Vector *BaseAdaptor_uncachedFetchAllByDbIDList(BaseAdaptor *ba, Vector *idList, Slice *slice);
SeqFeature *BaseAdaptor_uncachedFetchByDbID(BaseAdaptor *ba, IDType id);
Vector *BaseAdaptor_fetchAll(BaseAdaptor *ba);
int BaseAdaptor_hasNoIdCache(BaseAdaptor *ba);

#define BaseAdaptor_setSpeciesId(ba,val) (ba)->speciesId = (val)
#define BaseAdaptor_getSpeciesId(ba) (ba)->speciesId

#define BaseAdaptor_isMultiSpecies(ba) (0)
#endif
