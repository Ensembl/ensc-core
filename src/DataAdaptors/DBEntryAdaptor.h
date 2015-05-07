/*
 * Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

#ifndef __DBENTRYADAPTOR_H__
#define __DBENTRYADAPTOR_H__

#include "BaseAdaptor.h"
#include "AdaptorTypes.h"
#include "DBEntry.h"

struct DBEntryAdaptorStruct {
  BASEADAPTOR_DATA
};

IDType DBEntryAdaptor_exists(DBEntryAdaptor *dbea, DBEntry *dbe);
DBEntryAdaptor *DBEntryAdaptor_new(DBAdaptor *dba);
DBEntry *DBEntryAdaptor_fetchByDbID(DBEntryAdaptor *dbea, IDType dbID);
IDType DBEntryAdaptor_store(DBEntryAdaptor *dbea, DBEntry *exObj,
                         IDType ensObject, char *ensType, int ignoreRelease);
int DBEntryAdaptor_fetchAllByTranscript(DBEntryAdaptor *dbea, Transcript *trans);
int DBEntryAdaptor_fetchAllByGene(DBEntryAdaptor *dbea, Gene *gene);
Vector *DBEntryAdaptor_fetchAllByTranslation(DBEntryAdaptor *dbea, Translation *trans);



#endif
