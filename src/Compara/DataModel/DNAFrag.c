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

#define __DNAFRAG_MAIN__
#include "DNAFrag.h"
#undef __DNAFRAG_MAIN__
#include "StrUtil.h"
#include "BaseContig.h"
#include "DBAdaptor.h"
#include "RawContigAdaptor.h"
#include "SliceAdaptor.h"

DNAFrag *DNAFrag_new() {
  DNAFrag *df;

  if ((df = (DNAFrag *)calloc(1,sizeof(DNAFrag))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for df\n");
    return NULL;
  }

  df->objectType = CLASS_DNAFRAG;

  df->funcs = &dnaFragFuncs;

  Object_incRefCount(df);
  return df;
}

ECOSTRING DNAFrag_setName(DNAFrag *df, char *name) {
  EcoString_copyStr(ecoSTable,&(df->name),name,0);

  return df->name;
}

char *DNAFrag_setType(DNAFrag *df, char *type) {
  StrUtil_copyString(&(df->type),type,0);

  return df->type;
}

BaseContig *DNAFrag_getContig(DNAFrag *df) {

   if (!df->contig) {
     DBAdaptor *dba = GenomeDB_getDBAdaptor(DNAFrag_getGenomeDB(df));
     if (!strcmp(DNAFrag_getType(df),"RawContig")) {
       RawContigAdaptor *rca = DBAdaptor_getRawContigAdaptor(dba);
       df->contig = (BaseContig *)RawContigAdaptor_fetchByName(rca, DNAFrag_getName(df));
     } else if (!strcmp(DNAFrag_getType(df),"VirtualContig")) {
       fprintf(stderr, "VC type not implemented \n");
       //my ($chr,$start,$end) = split /\./, $self->name;
       //df->contig = $core_dbadaptor->get_SliceAdaptor->fetch_by_chr_start_end(chrName,start,end);
     } else if (!strcmp(DNAFrag_getType(df),"Chromosome")) {
       SliceAdaptor *sa = DBAdaptor_getSliceAdaptor(dba);
       df->contig = (BaseContig *)SliceAdaptor_fetchByChrName(sa, DNAFrag_getName(df));
     } else {
       fprintf(stderr, "Error: Can't fetch contig of %s with type %s\n",
               DNAFrag_getName(df), DNAFrag_getType(df));
     }
   }

   return df->contig;
}

void DNAFrag_free(DNAFrag *df) {
  Object_decRefCount(df);

  if (Object_getRefCount(df) > 0) {
    return;
  } else if (Object_getRefCount(df) < 0) {
    fprintf(stderr,"Error: Negative reference count for DNAFrag\n"
                   "       Freeing it anyway\n");
  }

  printf("BaseContig_free needs implementing in DNAFrag\n");
  //if (df->contig)   BaseContig_free(df->contig);
  if (df->genomeDB) GenomeDB_free(df->genomeDB);

  if (df->type) EcoString_freeStr(ecoSTable, df->type);
  if (df->name) EcoString_freeStr(ecoSTable, df->name);

  free(df);
}

