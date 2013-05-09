#ifndef __BASEFEATUREADAPTOR_H__
#define __BASEFEATUREADAPTOR_H__

#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "DBAdaptor.h"
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
  int startEqualsEnd; \
  long maxFeatureLen;
/*
  BaseFeatureAdaptor_StoreFunc store; \
  BaseFeatureAdaptor_GetTablesFunc getTables; \
  BaseFeatureAdaptor_GetColumnsFunc getColumns; \
  BaseFeatureAdaptor_ObjectsFromStatementHandleFunc objectsFromStatementHandle; \
  BaseFeatureAdaptor_FinalClauseFunc finalClause; \
  BaseFeatureAdaptor_LeftJoinFunc leftJoin; \
  BaseFeatureAdaptor_DefaultWhereClauseFunc defaultWhereClause;
*/

struct BaseFeatureAdaptorStruct {
  BASEFEATUREADAPTOR_DATA
};

#define BaseFeatureAdaptor_setStartEqualsEnd(bfa,val) (bfa)->startEqualsEnd = (val)
#define BaseFeatureAdaptor_getStartEqualsEnd(bfa) (bfa)->startEqualsEnd

#define BaseFeatureAdaptor_setMaxFeatureLength(bfa,val) (bfa)->maxFeatureLen = (val)
#define BaseFeatureAdaptor_getMaxFeatureLength(bfa) (bfa)->maxFeatureLen

#define BaseFeatureAdaptor_setSpeciesId(bfa,val) BaseAdaptor_setSpeciesId((bfa), (val))
#define BaseFeatureAdaptor_getSpeciesId(bfa) BaseAdaptor_getSpeciesId((bfa))

#define BaseFeatureAdaptor_fetchByDbID(bfa,val) BaseAdaptor_fetchByDbID((BaseFeatureAdaptor *)(bfa), (val))
/*
char *BaseFeatureAdaptor_finalClause(void);
char **BaseFeatureAdaptor_leftJoin(void);
char *BaseFeatureAdaptor_defaultWhereClause(void);
int BaseFeatureAdaptor_store(BaseFeatureAdaptor *bfa, Vector *features);
*/

void BaseFeatureAdaptor_init(BaseFeatureAdaptor *bfa, DBAdaptor *dba, int adaptorType);
void BaseFeatureAdaptor_clearSliceFeatureCache(BaseFeatureAdaptor *bfa);

Vector *BaseFeatureAdaptor_fetchAllBySlice(BaseFeatureAdaptor *bfa, Slice *slice, char *logicName);
Vector *BaseFeatureAdaptor_fetchAllBySliceAndScore(BaseFeatureAdaptor *bfa, Slice *slice,
                                                   double *scoreP, char *logicName);
Vector *BaseFeatureAdaptor_fetchAllBySliceConstraint(BaseFeatureAdaptor *bfa, Slice *slice, char *constraint, char *logicName);
Vector *BaseFeatureAdaptor_fetchAllByLogicName(BaseFeatureAdaptor *bfa, char *logicName);
Vector *BaseFeatureAdaptor_fetchAllByStableIdList(BaseFeatureAdaptor *bfa, Vector *ids, Slice *slice);
int BaseFeatureAdaptor_countBySliceConstraint(BaseFeatureAdaptor *bfa, Slice *slice, char *constraint, char *logicName);
Vector *BaseFeatureAdaptor_getAndFilterSliceProjections(BaseFeatureAdaptor *bfa, Slice *slice);
long *BaseFeatureAdaptor_generateFeatureBounds(BaseFeatureAdaptor *bfa, Slice *slice, int *nBound);
Vector *BaseFeatureAdaptor_getBySlice(BaseFeatureAdaptor *bfa, Slice *slice, char *origConstraint, char *queryType);
Vector *BaseFeatureAdaptor_sliceFetch(BaseFeatureAdaptor *bfa, Slice *slice, char *origConstraint);
Vector *BaseFeatureAdaptor_remap(BaseFeatureAdaptor *bfa, Vector *features, AssemblyMapper *mapper, Slice *slice);
char *BaseFeatureAdaptor_logicNameToConstraint(BaseFeatureAdaptor *bfa, char *constraint, char *logicName);
Vector *BaseFeatureAdaptor_listSeqRegionIds(BaseFeatureAdaptor *bfa, char *table);


/*
NameTableType *BaseFeatureAdaptor_getTables(void); 
char *BaseFeatureAdaptor_getColumns(void); 
Vector *BaseFeatureAdaptor_objectsFromStatementHandle(BaseFeatureAdaptor *bfa, StatementHandle *sth,AssemblyMapper *mapper, Slice *slice); 
*/





#endif
