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
  CLASS_TRANSCRIPT,
  CLASS_GENE,
  CLASS_FEATURESET,
  CLASS_NUMCLASS
} ClassType;

struct ClassHierarchyNodeStruct {
  int nSubClass;
  Class *class;
  ClassHierarchyNode **subClasses;
};

struct ClassStruct {
  ClassType type;
  char *name;
};


#endif
