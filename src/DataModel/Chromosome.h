#ifndef __CHROMOSOME_H__
#define __CHROMOSOME_H__

#include "DataModelTypes.h"
#include "Storable.h"
#include "Object.h"

struct ChromosomeStruct {
  OBJECT_DATA
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

ECOSTRING Chromosome_setName(Chromosome *c,char *name);
#define Chromosome_getName(c) (c)->name

void Chromosome_free(Chromosome *chr);

#endif
