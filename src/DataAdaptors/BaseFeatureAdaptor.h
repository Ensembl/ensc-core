#ifndef __BASEFEATUREADAPTOR_H__
#define __BASEFEATUREADAPTOR_H__

#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "Set.h"
#include "StatementHandle.h"

typedef int      (*BaseFeatureAdaptor_StoreFunc)(BaseFeatureAdaptor *bfa, Set *features);
typedef char *** (*BaseFeatureAdaptor_GetTablesFunc)(void);
typedef char *   (*BaseFeatureAdaptor_GetColumnsFunc)(void);
typedef Set *    (*BaseFeatureAdaptor_ObjectsFromStatementHandleFunc)(BaseFeatureAdaptor *bfa,StatementHandle *sth);

#define BASEFEATUREADAPTOR_DATA
  BASEADAPTOR_DATA \
  BaseFeatureAdaptor_StoreFunc store; \
  BaseFeatureAdaptor_GetTablesFunc getTables; \
  BaseFeatureAdaptor_GetColumnsFunc getColumns; \
  BaseFeatureAdaptor_ObjectsFromStatementHandleFunc objectsFromStatementHandle;

struct BaseFeatureAdaptorStruct {
  BASEFEATUREADAPTOR_DATA
};

#endif
