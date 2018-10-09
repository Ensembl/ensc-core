/*
 * Copyright [1999-2017] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

#ifndef __BASERODBTEST_H__
#define __BASERODBTEST_H__

#include "BaseTest.h"
#include "SliceAdaptor.h"

#include "DBAdaptor.h"

DBAdaptor *Test_initROEnsDB() {
  DBAdaptor *dba;

  dba = DBAdaptor_new("ensembldb.ensembl.org","anonymous",NULL,"homo_sapiens_core_86_38",3306,NULL);

  return dba;
}

Slice *Test_getStandardSlice(DBAdaptor *dba) {
  Slice *slice;
  SliceAdaptor *sa;

  sa = DBAdaptor_getSliceAdaptor(dba);

  slice = SliceAdaptor_fetchByRegion(sa,"chromosome","20",1000000,2600000,1,NULL,0);

  return slice;
}
#endif 
