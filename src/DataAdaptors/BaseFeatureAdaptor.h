#ifndef __BASEFEATUREADAPTOR_H__
#define __BASEFEATUREADAPTOR_H__

#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "Vector.h"
#include "RawContig.h"
#include "StatementHandle.h"
#include "Cache.h"
#include "Slice.h"
#include "AssemblyMapper.h"

typedef char * NameTableType[][2];

typedef int      (*BaseFeatureAdaptor_StoreFunc)(BaseFeatureAdaptor *bfa, Vector *features);
typedef NameTableType *(*BaseFeatureAdaptor_GetTablesFunc)(void);
typedef char *   (*BaseFeatureAdaptor_GetColumnsFunc)(void);
typedef Vector *    (*BaseFeatureAdaptor_ObjectsFromStatementHandleFunc)(BaseFeatureAdaptor *bfa,StatementHandle *sth,
                                                                      AssemblyMapper *mapper, Slice *slice);
typedef char *   (*BaseFeatureAdaptor_FinalClauseFunc)(void);
typedef char *   (*BaseFeatureAdaptor_DefaultWhereClauseFunc)(void);
typedef char **  (*BaseFeatureAdaptor_LeftJoinFunc)(void);


#define BASEFEATUREADAPTOR_DATA \
  BASEADAPTOR_DATA \
  Cache *sliceFeatureCache; \
  BaseFeatureAdaptor_StoreFunc store; \
  BaseFeatureAdaptor_GetTablesFunc getTables; \
  BaseFeatureAdaptor_GetColumnsFunc getColumns; \
  BaseFeatureAdaptor_ObjectsFromStatementHandleFunc objectsFromStatementHandle; \
  BaseFeatureAdaptor_FinalClauseFunc finalClause; \
  BaseFeatureAdaptor_LeftJoinFunc leftJoin; \
  BaseFeatureAdaptor_DefaultWhereClauseFunc defaultWhereClause;

struct BaseFeatureAdaptorStruct {
  BASEFEATUREADAPTOR_DATA
};

char *BaseFeatureAdaptor_finalClause(void);
char **BaseFeatureAdaptor_leftJoin(void);
char *BaseFeatureAdaptor_defaultWhereClause(void);
int BaseFeatureAdaptor_store(BaseFeatureAdaptor *bfa, Vector *features);

void BaseFeatureAdaptor_init(BaseFeatureAdaptor *bfa, DBAdaptor *dba, int adaptorType);
Vector *BaseFeatureAdaptor_genericFetch(BaseFeatureAdaptor *bfa, char *constraint,
                                     char *logicName, AssemblyMapper *mapper, Slice *slice);
SeqFeature *BaseFeatureAdaptor_fetchByDbID(BaseFeatureAdaptor *bfa, IDType dbID);
Vector *BaseFeatureAdaptor_fetchAllByRawContigConstraint(BaseFeatureAdaptor *bfa, RawContig *contig,
                                                      char *constraint, char *logicName);
Vector *BaseFeatureAdaptor_fetchAllByRawContig(BaseFeatureAdaptor *bfa, RawContig *contig,
                                            char *logicName);
Vector *BaseFeatureAdaptor_fetchAllByRawContigAndScore(BaseFeatureAdaptor *bfa, RawContig *contig,
                                                    double *scoreP, char *logicName);
Vector *BaseFeatureAdaptor_fetchAllBySlice(BaseFeatureAdaptor *bfa, Slice *slice,
                                        char *logicName);
Vector *BaseFeatureAdaptor_fetchAllBySliceAndScore(BaseFeatureAdaptor *bfa, Slice *slice,
                                                double *scoreP, char *logicName);
Vector *BaseFeatureAdaptor_fetchAllBySliceConstraint(BaseFeatureAdaptor *bfa, Slice *slice,
                                                  char *constraint, char *logicName);
int BaseFeatureAdaptor_remove(BaseFeatureAdaptor *bfa, SeqFeature *feature);
int BaseFeatureAdaptor_removeByRawContig(BaseFeatureAdaptor *bfa, RawContig *contig);
NameTableType *BaseFeatureAdaptor_getTables(void); 
char *BaseFeatureAdaptor_getColumns(void); 
Vector *BaseFeatureAdaptor_objectsFromStatementHandle(BaseFeatureAdaptor *bfa, StatementHandle *sth,AssemblyMapper *mapper, Slice *slice); 





#endif
