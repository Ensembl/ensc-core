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

#ifndef __GENOMEDB_H__
#define __GENOMEDB_H__

#include "ComparaDataModelTypes.h"
#include "EnsRoot.h"
#include "Storable.h"
#include "Vector.h"

OBJECTFUNC_TYPES(GenomeDB)

typedef struct GenomeDBFuncsStruct {
  OBJECTFUNCS_DATA(GenomeDB)
} GenomeDBFuncs;

#define FUNCSTRUCTTYPE GenomeDBFuncs
struct GenomeDBStruct {
  ENSROOT_DATA
  Storable st;
  char *name;
  DBAdaptor *dbAdaptor;
  char *assembly;
  IDType taxonId;
};
#undef FUNCSTRUCTTYPE

#define GenomeDB_setDbID(gd,id) Storable_setDbID(&((gd)->st),(id))
#define GenomeDB_getDbID(gd) Storable_getDbID(&((gd)->st))

#define GenomeDB_setAdaptor(gd,a) Storable_setAdaptor(&((gd)->st),(a))
#define GenomeDB_getAdaptor(gd) Storable_getAdaptor(&((gd)->st))

#define GenomeDB_setDBAdaptor(gd,a) (gd)->dbAdaptor = (a)
#define GenomeDB_getDBAdaptor(gd) (gd)->dbAdaptor

char *GenomeDB_setAssembly(GenomeDB *gd, char *assembly);
#define GenomeDB_getAssembly(gd) (gd)->assembly

char *GenomeDB_setName(GenomeDB *gd, char *name);
#define GenomeDB_getName(gd) (gd)->name

#define GenomeDB_setTaxonId(gd,t) (gd)->taxonId = (t)
#define GenomeDB_getTaxonId(gd) (gd)->taxonId

GenomeDB *GenomeDB_new();

Vector *GenomeDB_linkedGenomesByMethodLinkId(GenomeDB *gdb, IDType methodLinkId);

void GenomeDB_free(GenomeDB *gd);

#ifdef __GENOMEDB_MAIN__
  GenomeDBFuncs
    genomeDBFuncs = {
                     GenomeDB_free
                    };
#else
  extern GenomeDBFuncs enomeDBFuncs;
#endif



#endif
