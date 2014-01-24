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

#include "SyntenyRegionAdaptor.h"

SyntenyRegionAdaptor *SyntenyRegionAdaptor_new(ComparaDBAdaptor *dba) {
  SyntenyRegionAdaptor *sra;

  if ((sra = (SyntenyRegionAdaptor *)calloc(1,sizeof(SyntenyRegionAdaptor))) == NULL) {
    fprintf(stderr,"Error: Failed allocating sra\n");
    exit(1);
  }
  BaseComparaAdaptor_init((BaseComparaAdaptor *)sra, dba, SYNTENYREGION_ADAPTOR);



  return sra;
}

SyntenyRegion *SyntenyRegionAdaptor_fetchByDbID(SyntenyRegionAdaptor *sra, IDType dbID) {
  StatementHandle *sth;
  ResultRow *row;
  char qStr[1024];
  SyntenyRegion *sr;
  
  if (!dbID) {
    fprintf(stderr,"Error: SyntenyRegionAdaptor_fetchByDbID with no dbID!\n");
    exit(1);
  }

  sprintf(qStr,
           "select synteny_cluster_id,"
           "       dnafrag_id,"
           "       seq_start,"
           "       seq_end"
           " from synteny_region "
           " where synteny_region_id = " IDFMTSTR,
             dbID);
           
  sth = sra->prepare((BaseAdaptor *)sra, qStr, strlen(qStr));
  sth->execute(sth);

  if (!(row = sth->fetchRow(sth))) {
    fprintf(stderr,"Error: No such dbID " IDFMTSTR "\n",dbID);
    exit(1);
  }

  sr =  SyntenyRegionAdaptor_newRegionFromArray(sra, 
                                        dbID,
                                        row->getLongLongAt(row,0),
                                        row->getLongLongAt(row,1), 
                                        row->getIntAt(row,2),
                                        row->getIntAt(row,3));
  sth->finish(sth);

  return sr;
}

Vector *SyntenyRegionAdaptor_fetchByClusterId(SyntenyRegionAdaptor *sra, IDType clusterId) {
  char qStr[256];
  StatementHandle *sth;
  ResultRow *row;
  Vector *out;

  if (!clusterId) {
    fprintf(stderr, "Error: fetch_by_cluster_id with no cluster_id!\n");
    exit(1);
  }

  sprintf(qStr,
            "select synteny_region_id,"
            "       dnafrag_id,"
            "       seq_start,"
            "       seq_end "
            " from synteny_region "
            " where synteny_cluster_id = " IDFMTSTR, clusterId);

  sth = sra->prepare((BaseAdaptor *)sra, qStr, strlen(qStr));
  sth->execute(sth);

  out = Vector_new();
  while ((row = sth->fetchRow(sth))) {
    SyntenyRegion *sr = SyntenyRegionAdaptor_newRegionFromArray(sra, 
                              row->getLongLongAt(row,0),
                              clusterId,
                              row->getLongLongAt(row,1), 
                              row->getIntAt(row,2),
                              row->getIntAt(row,3));
    
    Vector_addElement(out,sr);
  }
  sth->finish(sth);

  return out;
}

SyntenyRegion *SyntenyRegionAdaptor_newRegionFromArray(SyntenyRegionAdaptor *sra, IDType dbID, 
                                            IDType cluster, IDType dnaFrag,int start,int end) {
  SyntenyRegion *region;

  region = SyntenyRegion_new();

  SyntenyRegion_setClusterId(region, cluster);
  SyntenyRegion_setDNAFragId(region, dnaFrag);
  SyntenyRegion_setStart(region, start);
  SyntenyRegion_setEnd(region, end);
  SyntenyRegion_setAdaptor(region, sra);
  SyntenyRegion_setDbID(region, dbID);

  return region;
}

IDType SyntenyRegionAdaptor_store(SyntenyRegionAdaptor *sra, IDType clusterId, SyntenyRegion *region) {
  IDType regionId;
  StatementHandle *sth;
  char qStr[1024];
  
  if (!region) {
    fprintf(stderr, "Error: region is not a SyntenyRegion\n");
    exit(1);
  }

  sprintf(qStr, "insert into synteny_region (synteny_cluster_id,dnafrag_id,seq_start,seq_end)"
                " VALUES (" IDFMTSTR ", " IDFMTSTR ", %d, %d)", 
                clusterId,
                SyntenyRegion_getDNAFragId(region),
                SyntenyRegion_getSeqStart(region),
                SyntenyRegion_getSeqEnd(region));
  sth = sra->prepare((BaseAdaptor *)sra, qStr, strlen(qStr));
   
  sth->execute(sth);

  regionId = sth->getInsertId(sth);
   
  SyntenyRegion_setDbID(region, regionId);
  SyntenyRegion_setAdaptor(region, sra);
   
  return regionId;
}
