#include "Class.h"
#include "EnsC.h"

#include <stdlib.h>
#include <string.h>

/* !!! Three !!! steps to adding an class */
/* Step 1: Add to hierarchyString */
/* Step 2: Add to classArray */
/* Step 3: Add to ClassType in Class.h */

/* Step 1 here */
char * hierarchyString =  
  " VECTOR\n"
  " STATEMENTHANDLE\n"
  "  MYSQLSTATEMENTHANDLE\n"
  " RESULTROW\n"
  "  MYSQLRESULTROW\n"
  " ENSROOT\n"
  "  DBENTRY\n"
  "  ANALYSIS\n"
  "  CLONE\n"
  "  CHROMOSOME\n"
  "  DNAFRAG\n"
  "  SYNTENYREGION\n"
  "  GENOMICALIGN\n"
  "  HOMOLOGY\n"
  "  GENOMEDB\n"
  "  SEQUENCE\n"
  "   BASECONTIG\n"
  "    RAWCONTIG\n"
  "    SLICE\n"
  "  REPEATCONSENSUS\n"
  "  TRANSLATION\n"
  "  SEQFEATURE\n"
  "   ANNOTATEDSEQFEATURE\n"
  "    EXON\n"
  "     STICKYEXON\n"
  "    TRANSCRIPT\n"
  "    GENE\n"
  "   PREDICTIONTRANSCRIPT\n"
  "   SIMPLEFEATURE\n"
  "   REPEATFEATURE\n"
  "   FEATUREPAIR\n"
  "    BASEALIGNFEATURE\n"
  "     DNADNAALIGNFEATURE\n"
  "     DNAPEPALIGNFEATURE\n"
  "  SPECIES\n"
  "  BASEASSEMBLYMAPPER\n"
  "   ASSEMBLYMAPPER\n"
  "   CHAINEDASSEMBLYMAPPER\n"
  "   TOPLEVELASSEMBLYMAPPER\n"
  "  COORDSYSTEM\n"
  "  ATTRIBUTE\n"
  "  SEQEDIT\n";

/* Step 2 here */
Class classArray[CLASS_NUMCLASS] = 
 {{CLASS_NONE, "NONE"},
  {CLASS_OBJECT, "OBJECT"},
  {CLASS_STATEMENTHANDLE, "STATEMENTHANDLE"},
  {CLASS_MYSQLSTATEMENTHANDLE, "MYSQLSTATEMENTHANDLE"},
  {CLASS_RESULTROW, "RESULTROW"},
  {CLASS_MYSQLRESULTROW, "MYSQLRESULTROW"},
  {CLASS_SEQFEATURE, "SEQFEATURE"},
  {CLASS_SIMPLEFEATURE, "SIMPLEFEATURE"},
  {CLASS_REPEATFEATURE, "REPEATFEATURE"},
  {CLASS_FEATUREPAIR, "FEATUREPAIR"},
  {CLASS_BASEALIGNFEATURE, "BASEALIGNFEATURE"},
  {CLASS_DNADNAALIGNFEATURE, "DNADNAALIGNFEATURE"},
  {CLASS_DNAPEPALIGNFEATURE, "DNAPEPALIGNFEATURE"},
  {CLASS_EXON, "EXON"},
  {CLASS_STICKYEXON, "STICKYEXON"},
  {CLASS_GENE, "GENE"},
  {CLASS_TRANSCRIPT, "TRANSCRIPT"},
  {CLASS_TRANSLATION, "TRANSLATION"},
  {CLASS_REPEATCONSENSUS, "REPEATCONSENSUS"},
  {CLASS_ENSROOT, "ENSROOT"},
  {CLASS_DBENTRY, "DBENTRY"},
  {CLASS_ANALYSIS, "ANALYSIS"},
  {CLASS_RAWCONTIG, "RAWCONTIG"},
  {CLASS_BASECONTIG, "BASECONTIG"},
  {CLASS_SEQUENCE, "SEQUENCE"},
  {CLASS_SLICE, "SLICE"},
  {CLASS_CHROMOSOME, "CHROMOSOME"},
  {CLASS_CLONE, "CLONE"},
  {CLASS_DNAFRAG, "DNAFRAG"},
  {CLASS_GENOMEDB, "GENOMEDB"},
  {CLASS_GENOMICALIGN, "GENOMICALIGN"},
  {CLASS_HOMOLOGY, "HOMOLOGY"},
  {CLASS_SYNTENYREGION, "SYNTENYREGION"},
  {CLASS_ANNOTATEDSEQFEATURE, "ANNOTATEDSEQFEATURE"},
  {CLASS_PREDICTIONTRANSCRIPT, "PREDICTIONTRANSCRIPT"},
  {CLASS_FEATURESET, "FEATURESET"},
  {CLASS_VECTOR, "VECTOR"},
  {CLASS_SPECIES, "SPECIES"},
  {CLASS_BASEASSEMBLYMAPPER, "BASEASSEMBLYMAPPER"},
  {CLASS_ASSEMBLYMAPPER, "ASSEMBLYMAPPER"},
  {CLASS_CHAINEDASSEMBLYMAPPER, "CHAINEDASSEMBLYMAPPER"},
  {CLASS_TOPLEVELASSEMBLYMAPPER, "TOPLEVELASSEMBLYMAPPER"},
  {CLASS_COORDSYSTEM, "COORDSYSTEM"},
  {CLASS_ATTRIBUTE, "ATTRIBUTE"},
  {CLASS_SEQEDIT, "SEQEDIT"}
  };

ClassHierarchyNode *root = NULL;

void Class_initHierarchy(ClassHierarchyNode *ch,  char **chPP,int *depth);
Class *Class_findByName(char *name);
int Class_searchHierarchyForParentChildPair(ClassHierarchyNode *chn, ClassType parentType, 
                                             ClassType descType, int belowParent, int *depth);

ClassHierarchyNode *ClassHierarchyNode_new(char *className) {
  ClassHierarchyNode *chn;
  if ((chn = (ClassHierarchyNode *) calloc(1,sizeof(ClassHierarchyNode))) == NULL) {
    fprintf(stderr, "Error: Failed allocating ClassHierarchyNode\n");
    exit(1);
  }
  chn->class = Class_findByName(className);
  return chn;
}

ClassHierarchyNode *ClassHierarchyNode_addSubClass(ClassHierarchyNode *parent,ClassHierarchyNode *child) {
  parent->nSubClass++;
  if ((parent->subClasses = (ClassHierarchyNode **)realloc(parent->subClasses,parent->nSubClass*sizeof(ClassHierarchyNode))) == NULL) {
    fprintf(stderr, "Error: Failed allocating ClassHierarchyNode subClasses\n");
    exit(1);
  }
  parent->subClasses[parent->nSubClass-1] = child;
  return child;
}

int Class_assertType(ClassType wantedType, ClassType actualType) {
  if (!Class_isDescendent(wantedType,actualType)) {
    fprintf(stderr,"Class type assertion failed for %d and %d\n",wantedType,actualType);
    exit(1);
  }
  return 1;
}

int Class_isDescendent(ClassType parentType, ClassType descType) {
  int sDepth = 0;

  if (root == NULL) {
    int ihDepth = 1;
    char *chP = hierarchyString;

    root = ClassHierarchyNode_new("OBJECT");
    Class_initHierarchy(root, &chP, &ihDepth);
  }

  if (descType == parentType) {
    return TRUE;
  } 

  return Class_searchHierarchyForParentChildPair(root, parentType, descType, FALSE, &sDepth);
}

int Class_searchHierarchyForParentChildPair(ClassHierarchyNode *chn, ClassType parentType, 
                                             ClassType descType, int belowParent, int *depth) {
  int i;
  int isParent = 0;

  (*depth)++;

  if (chn->class->type == parentType) {
    if (belowParent) {
      fprintf(stderr,"ERROR: Class type %s is parent to itself\n", chn->class->name);
      exit(1);
    }
    belowParent = TRUE;
    isParent = 1;
  } else if (chn->class->type == descType && belowParent == TRUE) {
    (*depth)--;
    return 1;
  }

  for (i=0;i<chn->nSubClass;i++) {
    if (Class_searchHierarchyForParentChildPair(chn->subClasses[i], parentType, descType, belowParent, depth)) {
      (*depth)--;
      return 1;
    }
  }

  if (isParent) {
    belowParent = FALSE;
  }
 
  (*depth)--;

  return 0;
}

#define MAXCLASSTYPE 1024
void Class_initHierarchy(ClassHierarchyNode *chn, char **chPP,int *depth) {
  Class **addedClasses[MAXCLASSTYPE];
  int nAdded = 0;
  char clsName[1024];
  int len=0;
  int curDepth;
  int i;
  char *chP = *chPP;

  while (*chP != '\0') {
    //printf("%c",*chP);
    if (*chP == '\n') {
      clsName[len] = '\0';
      chP++;

      curDepth = 0;
      /* find depth */
      while (clsName[curDepth] == ' ') {
        curDepth++;
      }

      //printf("depth = %d curDepth = %d\n",*depth, curDepth);
      if (*depth == curDepth) { /* If its same - add */

        *chPP = chP;
        ClassHierarchyNode *newNode = ClassHierarchyNode_new(&clsName[curDepth]);

        ClassHierarchyNode_addSubClass(chn, newNode);
        //printf(" same depth %s\n",clsName);

      } else if (*depth < curDepth) { /* If its greater - start new sub classes (recurse) */

        (*depth)++;
        //printf("descending...\n");
        Class_initHierarchy(chn->subClasses[chn->nSubClass-1], chPP, depth);
        chP = *chPP;
      } else { /* If its less - return */
        (*depth)--;
        //printf("returning ...\n");
        return;
      }
      len=0;
    } else {
      clsName[len++] = *chP;
      chP++;
    }
  }
}

int Class_checkIfAdded(Class **cls, int nObj, Class *newOne) {
  int i;

  for (i=0;i<nObj;i++) {
    if (cls[i] == newOne) {
      return 1;
    }
  }
  return 0;
}

/* searches the global list of class types */
Class *Class_findByName(char *name) {
  int i;

  for (i=0;i<CLASS_NUMCLASS;i++) {
    if (!strcmp(classArray[i].name, name)) {
      return &classArray[i];
    }
  } 
  fprintf(stderr,"ERROR: Class type %s not found in class types array\n",name);
  exit(1);
}

Class *Class_findByType(ClassType type) {
  int i;

  for (i=0;i<CLASS_NUMCLASS;i++) {
    if (type == classArray[i].type) {
      return &classArray[i];
    }
  } 
  fprintf(stderr,"ERROR: Class type %d not found in class types array\n",type);
  exit(1);
}

