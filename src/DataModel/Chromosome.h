#ifndef __CHROMOSOME_H__
#define __CHROMOSOME_H__

#include "DataModelTypes.h"
#include "Storable.h"

struct ChromosomeStruct {
  Storable st;
  int length;
  char *name;
};

Chromosome *Chromosome_new(void);

#define Chromosome_setDbID(c,dbID) Storable_setDbID(&((c)->st),(dbID))
#define Chromosome_getDbID(c) Storable_getDbID(&((c)->st))

#define Chromosome_setAdaptor(c,ad) Storable_setAdaptor(&((c)->st),(ad))
#define Chromosome_getAdaptor(c) Storable_getAdaptor(&((c)->st))

#define Chromosome_setLength(c,len) (c)->length = (len)
#define Chromosome_getLength(c) (c)->length

#define Chromosome_setName(c,n) (c)->name = (n)
#define Chromosome_getName(c) (c)->name

#endif
