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

#ifndef __ANALYSISADAPTOR_H__
#define __ANALYSISADAPTOR_H__

#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "Analysis.h"
#include "IDHash.h"
#include "StringHash.h"

struct AnalysisAdaptorStruct {
  BASEADAPTOR_DATA
  IDHash *analCache;
  StringHash *logicNameCache;
};

IDType AnalysisAdaptor_analysisExists(AnalysisAdaptor *aa, Analysis *anal);
AnalysisAdaptor *AnalysisAdaptor_new(DBAdaptor *dba);
Analysis *AnalysisAdaptor_fetchByDbID(AnalysisAdaptor *aa, IDType dbID);
Analysis *AnalysisAdaptor_fetchByLogicName(AnalysisAdaptor *aa, char *logicName);
Analysis *AnalysisAdaptor_analysisFromRow(AnalysisAdaptor *aa, ResultRow *row);
Analysis **AnalysisAdaptor_fetchAll(AnalysisAdaptor *aa);
IDType AnalysisAdaptor_store(AnalysisAdaptor *aa, Analysis *analysis);






#endif
