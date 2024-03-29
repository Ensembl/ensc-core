/*
 * Copyright [1999-2024] EMBL-European Bioinformatics Institute
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

#ifndef __BASECOMPARADBTEST_H__
#define __BASECOMPARADBTEST_H__

#include "BaseTest.h"

#include "ComparaDBAdaptor.h"
#include "DBAdaptor.h"

ComparaDBAdaptor *Test_initComparaDB() {
  ComparaDBAdaptor *cdba;

  cdba = ComparaDBAdaptor_new("ensembldb.ensembl.org","anonymous",NULL,"ensembl_compara_24_1",4306,"Tests/Compara.conf");

  return cdba;
}

Slice *Test_getStandardSlice(ComparaDBAdaptor *cdba) {
  Slice *slice;
  SliceAdaptor *sa;
  DBAdaptor *dba = ComparaDBAdaptor_getDBAdaptor(cdba,"Homo sapiens", "NCBI34");

  sa = DBAdaptor_getSliceAdaptor(dba);

  slice = SliceAdaptor_fetchByRegion(sa,"chromosome","1",500000,1000000,1,NULL,0);

  return slice;
}
#endif 
