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

#ifndef __BASERODBTEST_H__
#define __BASERODBTEST_H__

#include "BaseTest.h"
#include "SliceAdaptor.h"

#include "DBAdaptor.h"

DBAdaptor *Test_initROEnsDB() {
  DBAdaptor *dba;

//  dba = DBAdaptor_new("kaka.sanger.ac.uk","anonymous",NULL,"homo_sapiens_core_14_31",3306,NULL);
//  dba = DBAdaptor_new("ecs2b.internal.sanger.ac.uk","ensro",NULL,"human_NCBI33_raw",3306,NULL);
//  dba = DBAdaptor_new("127.0.0.1","ensro",NULL,"human_NCBI33_raw",13302,NULL);
//  dba = DBAdaptor_new("127.0.0.1","root",NULL,"test_core",3306,NULL);
//  dba = DBAdaptor_new("127.0.0.1","root",NULL,"homo_sapiens_core_18_34",3306,NULL);
//  dba = DBAdaptor_new("127.0.0.1","root",NULL,"steve_feb04_comp",3306,NULL);
//  dba = DBAdaptor_new("ens-livemirror.internal.sanger.ac.uk","ensadmin","ensembl","homo_sapiens_core_71_37",3306,NULL);
  dba = DBAdaptor_new("127.0.0.1","ensadmin","ensembl","homo_sapiens_core_72_37",13390,NULL);
//  dba = DBAdaptor_new("genebuild2.internal.sanger.ac.uk","ensadmin","ensembl","steve_hs_testdb",3306,NULL);
//  dba = DBAdaptor_new("127.0.0.1","ensadmin","ensembl","steve_hs_testdb",3306,NULL);

  return dba;
}

Slice *Test_getStandardSlice(DBAdaptor *dba) {
  Slice *slice;
  SliceAdaptor *sa;

  sa = DBAdaptor_getSliceAdaptor(dba);

//  slice = SliceAdaptor_fetchByRegion(sa,"chromosome","2",10000000,100000000,1,NULL,0);
//  slice = SliceAdaptor_fetchByRegion(sa,"chromosome","1",1000000,260000000,1,NULL,0);
//  slice = SliceAdaptor_fetchByRegion(sa,"chromosome","2",1000000,260000000,1,NULL,0);
  //PAR + boundary cross
  //slice = SliceAdaptor_fetchByRegion(sa,"chromosome","Y",1000000,4000000,1,NULL,0);
//  slice = SliceAdaptor_fetchByRegion(sa,"chromosome","1",1000000,2000000,1,NULL,0);
  slice = SliceAdaptor_fetchByRegion(sa,"chromosome","20",1000000,260000000,1,NULL,0);

  return slice;
}
#endif 
