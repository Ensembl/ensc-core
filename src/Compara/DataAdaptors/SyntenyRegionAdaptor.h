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

#ifndef __SYNTENYREGIONADAPTOR_H__
#define __SYNTENYREGIONADAPTOR_H__

#include "BaseComparaAdaptor.h"
#include "ComparaAdaptorTypes.h"
#include "SyntenyRegion.h"
#include "Vector.h"


struct SyntenyRegionAdaptorStruct {
  BASECOMPARAADAPTOR_DATA
};

SyntenyRegionAdaptor *SyntenyRegionAdaptor_new(ComparaDBAdaptor *dba);
SyntenyRegion *SyntenyRegionAdaptor_newRegionFromArray(SyntenyRegionAdaptor *sra, IDType dbID,
                                            IDType cluster, IDType dnaFrag,int start,int end);
SyntenyRegion *SyntenyRegionAdaptor_fetchByDbID(SyntenyRegionAdaptor *sra, IDType dbID);

Vector *SyntenyRegionAdaptor_fetchByClusterId(SyntenyRegionAdaptor *sra, IDType clusterId);
IDType SyntenyRegionAdaptor_store(SyntenyRegionAdaptor *sra, IDType clusterId, SyntenyRegion *region);



#endif
