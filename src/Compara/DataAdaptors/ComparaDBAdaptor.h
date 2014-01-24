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

#ifndef __COMPARADBADAPTOR_H__
#define __COMPARADBADAPTOR_H__

#include "BaseDBAdaptor.h"
#include "DBConnection.h"
#include "ComparaAdaptorTypes.h"
#include "EnsC.h"
#include "StringHash.h"

struct ComparaDBAdaptorStruct {
  BASEDBADAPTOR_DATA
  StringHash *genomes;
};

ComparaDBAdaptor *ComparaDBAdaptor_new(char *host, char *user, char *pass, char *dbname,
                                       unsigned int port, char *confFile);

DBAdaptor *ComparaDBAdaptor_getDBAdaptor(ComparaDBAdaptor *cdba, char *species, char *assembly);
void ComparaDBAdaptor_addDBAdaptor(ComparaDBAdaptor *cdba, DBAdaptor *dba);
DBAdaptor *ComparaDBAdaptor_getDBAdaptor(ComparaDBAdaptor *cdba, char *species, char *assembly);
SyntenyAdaptor *ComparaDBAdaptor_getSyntenyAdaptor(ComparaDBAdaptor *cdba);
GenomeDBAdaptor *ComparaDBAdaptor_getGenomeDBAdaptor(ComparaDBAdaptor *cdba);
DNAFragAdaptor *ComparaDBAdaptor_getDNAFragAdaptor(ComparaDBAdaptor *cdba);
GenomicAlignAdaptor *ComparaDBAdaptor_getGenomicAlignAdaptor(ComparaDBAdaptor *cdba);
HomologyAdaptor *ComparaDBAdaptor_getHomologyAdaptor(ComparaDBAdaptor *cdba);
ComparaDNAAlignFeatureAdaptor *ComparaDBAdaptor_getComparaDNAAlignFeatureAdaptor(ComparaDBAdaptor *cdba);
MetaContainer *ComparaDBAdaptor_getMetaContainer(ComparaDBAdaptor *cdba);




#endif
