#ifndef __ANALYSIS_H__
#define __ANALYSIS_H__

#include "Storable.h"

typedef struct AnalysisStruct Analysis;

struct AnalysisStruct {
  Storable st;
  char *db;
  int   dbVersion;
  char *dbFile;
  char *program;
  int   programVersion;
  char *programFile;
  char *gffSource;
  char *gffFeature;
  char *module;
  int   moduleVersion;
  char *parameters;
  char *created;
  char *logicName;
};

Analysis *Analysis_new(void);

#define Analysis_setDbID(a,dbID) Storable_setDbID(&((a)->st),dbID)
#define Analysis_getDbID(a) Storable_getDbID(&((a)->st))

#define Analysis_setAdaptor(a,ad) Storable_setAdaptor(&((a)->st),ad)

#define Analysis_setLogicName(a,lname) (a)->logicName = (lname)
#define Analysis_getLogicName(a) (a)->logicName

#define Analysis_setProgram(a,prog) a->program = prog
#define Analysis_setProgramVersion(a,pver) a->programVersion = pver
#define Analysis_setProgramFile(a,pfile) a->programFile = pfile
#define Analysis_setDb(a,d) a->db = d
#define Analysis_setDbVersion(a,dver) a->dbVersion = dver
#define Analysis_setDbFile(a,dfile) a->dbFile = dfile
#define Analysis_setModule(a,mod) a->module = mod
#define Analysis_setModuleVersion(a,mver) a->moduleVersion = mver
#define Analysis_setGFFSource(a,gsrc) a->gffSource = gsrc
#define Analysis_setGFFFeature(a,gfeat) a->gffFeature = gfeat
#define Analysis_setCreated(a,c) a->created = c
#define Analysis_setParameters(a,params) a->parameters = params


#endif
