#ifndef __BASEFEATUREADAPTOR_H__
#define __BASEFEATUREADAPTOR_H__

#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "Set.h"
#include "StatementHandle.h"
#include "Cache.h"

typedef int      (*BaseFeatureAdaptor_StoreFunc)(BaseFeatureAdaptor *bfa, Set *features);
typedef char *** (*BaseFeatureAdaptor_GetTablesFunc)(void);
typedef char *   (*BaseFeatureAdaptor_GetColumnsFunc)(void);
typedef Set *    (*BaseFeatureAdaptor_ObjectsFromStatementHandleFunc)(BaseFeatureAdaptor *bfa,StatementHandle *sth);
typedef char *   (*BaseFeatureAdaptor_FinalClauseFunc)(void);
typedef char *   (*BaseFeatureAdaptor_DefaultWhereClauseFunc)(void);
typedef char **  (*BaseFeatureAdaptor_LeftJoinFunc)(void);


#define BASEFEATUREADAPTOR_DATA
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
int BaseFeatureAdaptor_store(BaseFeatureAdaptor *bfa, Set *features);

#endif
