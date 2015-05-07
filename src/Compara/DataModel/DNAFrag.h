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

#ifndef __DNAFRAG_H__
#define __DNAFRAG_H__

#include "ComparaDataModelTypes.h"
#include "EnsRoot.h"
#include "BaseContig.h"
#include "Storable.h"
#include "GenomeDB.h"


OBJECTFUNC_TYPES(DNAFrag)

typedef struct DNAFragFuncsStruct {
  OBJECTFUNCS_DATA(DNAFrag)
} DNAFragFuncs;

#define FUNCSTRUCTTYPE DNAFragFuncs
struct DNAFragStruct {
  ENSROOT_DATA
  Storable st;
  ECOSTRING name;
  BaseContig *contig;
  int start;
  int end;
  GenomeDB *genomeDB;
  ECOSTRING type;
};
#undef FUNCSTRUCTTYPE

BaseContig *DNAFrag_getContig(DNAFrag *df);

#define DNAFrag_setStart(df,s) (df)->start = (s)
#define DNAFrag_getStart(df) (df)->start

#define DNAFrag_setEnd(df,e) (df)->end = (e)
#define DNAFrag_getEnd(df) (df)->end

char *DNAFrag_setType(DNAFrag *df, char *type);
#define DNAFrag_getType(df) (df)->type

#define DNAFrag_setGenomeDB(df,g) (df)->genomeDB = (g)
#define DNAFrag_getGenomeDB(df) (df)->genomeDB

#define DNAFrag_setDbID(df,id) Storable_setDbID(&((df)->st),(id))
#define DNAFrag_getDbID(df) Storable_getDbID(&((df)->st))

#define DNAFrag_setAdaptor(df,a) Storable_setAdaptor(&((df)->st),(a))
#define DNAFrag_getAdaptor(df) Storable_getAdaptor(&((df)->st))

char *DNAFrag_setName(DNAFrag *df, char *name);
#define DNAFrag_getName(df) (df)->name

DNAFrag *DNAFrag_new();

void DNAFrag_free(DNAFrag *df);

#ifdef __DNAFRAG_MAIN__
  DNAFragFuncs
    dnaFragFuncs = {
                    DNAFrag_free
                   };
#else
  extern DNAFragFuncs dnaFragFuncs;
#endif


#endif
