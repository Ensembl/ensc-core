/*
 * Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
 * Copyright [2016-2017] EMBL-European Bioinformatics Institute
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

#ifndef __ANALYSIS_H__
#define __ANALYSIS_H__

#include "DataModelTypes.h"
#include "Storable.h"

#include "EnsRoot.h"

OBJECTFUNC_TYPES(Analysis)

typedef struct AnalysisFuncsStruct {
  OBJECTFUNCS_DATA(Analysis)
} AnalysisFuncs;

#define FUNCSTRUCTTYPE AnalysisFuncs
struct AnalysisStruct {
  ENSROOT_DATA
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
#undef FUNCSTRUCTTYPE

Analysis *Analysis_new(void);

#define Analysis_isStored(a, db) Storable_isStored(&((a)->st), (db))

#define Analysis_setDbID(a,dbID) Storable_setDbID(&((a)->st),dbID)
#define Analysis_getDbID(a) Storable_getDbID(&((a)->st))

#define Analysis_setAdaptor(a,ad) Storable_setAdaptor(&((a)->st),ad)
#define Analysis_getAdaptor(a) Storable_getAdaptor(&((a)->st))

#define Analysis_setLogicName(a,lname) (a)->logicName = (lname)
#define Analysis_getLogicName(a) (a)->logicName

#define Analysis_setProgram(a,prog) a->program = prog
#define Analysis_getProgram(a) a->program

#define Analysis_setProgramVersion(a,pver) a->programVersion = pver
#define Analysis_getProgramVersion(a) a->programVersion

#define Analysis_setProgramFile(a,pfile) a->programFile = pfile
#define Analysis_getProgramFile(a) a->programFile

#define Analysis_setDb(a,d) a->db = d
#define Analysis_getDb(a) a->db

#define Analysis_setDbVersion(a,dver) a->dbVersion = dver
#define Analysis_getDbVersion(a) a->dbVersion

#define Analysis_setDbFile(a,dfile) a->dbFile = dfile
#define Analysis_getDbFile(a) a->dbFile

#define Analysis_setModule(a,mod) a->module = mod
#define Analysis_getModule(a) a->module

#define Analysis_setModuleVersion(a,mver) a->moduleVersion = mver
#define Analysis_getModuleVersion(a) a->moduleVersion

#define Analysis_setGFFSource(a,gsrc) a->gffSource = gsrc
#define Analysis_getGFFSource(a) a->gffSource

#define Analysis_setGFFFeature(a,gfeat) a->gffFeature = gfeat
#define Analysis_getGFFFeature(a) a->gffFeature

#define Analysis_setCreated(a,c) a->created = c
#define Analysis_getCreated(a) a->created

#define Analysis_setParameters(a,params) a->parameters = params
#define Analysis_getParameters(a) a->parameters

int Analysis_compare(Analysis *a, Analysis *b);

void Analysis_free(Analysis *anal);

#ifdef __ANALYSIS_MAIN__
  AnalysisFuncs
    analysisFuncs = {
                    Analysis_free,
                    NULL, // shallowCopy
                    NULL  // deepCopy
                   };
#else
  extern AnalysisFuncs analysisFuncs;
#endif



#endif
