#ifndef __DNAFRAG_H__
#define __DNAFRAG_H__

#include "ComparaDataModelTypes.h"
#include "EnsRoot.h"

struct DNAFragStruct {
  ENSROOT_DATA
  Storable st;
  char *name;
  BaseContig *contig;
  int start;
  int end;
  char *genomeDB;
  char *type;
};

#endif
