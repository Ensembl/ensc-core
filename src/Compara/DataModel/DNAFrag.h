#ifndef __DNAFRAG_H__
#define __DNAFRAG_H__

#include "ComparaDataModelTypes.h"
#include "EnsRoot.h"
#include "BaseContig.h"
#include "Storable.h"
#include "GenomeDB.h"

#define FUNCSTRUCTTYPE NoTypeFuncs
struct DNAFragStruct {
  ENSROOT_DATA
  Storable st;
  char *name;
  BaseContig *contig;
  int start;
  int end;
  GenomeDB *genomeDB;
  char *type;
};
#undef FUNCSTRUCTTYPE

#define DNAFrag_setContig(df,c) (df)->contig = (c)
#define DNAFrag_getContig(df) (df)->contig

#define DNAFrag_setStart(df,s) (df)->start = (s)
#define DNAFrag_getStart(df) (df)->start

#define DNAFrag_setEnd(df,e) (df)->end = (e)
#define DNAFrag_getEnd(df) (df)->end

#define DNAFrag_setType(df,t) (df)->type = (t)
#define DNAFrag_getType(df) (df)->type

#define DNAFrag_setGenomeDB(df,g) (df)->genomeDB = (g)
#define DNAFrag_getGenomeDB(df) (df)->genomeDB

#define DNAFrag_setDbID(df,id) Storable_setDbID(&((df)->st),(id))
#define DNAFrag_getDbID(df) Storable_getDbID(&((df)->st))

#define DNAFrag_setAdaptor(df,a) Storable_setAdaptor(&((df)->st),(a))
#define DNAFrag_getAdaptor(df) Storable_getAdaptor(&((df)->st))

DNAFrag *DNAFrag_new();

#endif
