#ifndef __GENOMICRANGE_H__
#define __GENOMICRANGE_H__

#include "DataModelTypes.h"

struct GenomicRangeStruct {
  char *chrName;
  int   chrStart;
  int   chrEnd;
  long  chrId;
};

GenomicRange *GenomicRange_new();

#define GenomicRange_setChrStart(gr,s) (gr)->chrStart = (s)
#define GenomicRange_getChrStart(gr) (gr)->chrStart

#define GenomicRange_setChrEnd(gr,e) (gr)->chrEnd = (e)
#define GenomicRange_getChrEnd(gr) (gr)->chrEnd

char *setChrName(GenomicRange *gr, char *chrName);
#define GenomicRange_getChrName(gr) (gr)->chrName

#define GenomicRange_setChrId(gr,i) (gr)->chrId = (i)
#define GenomicRange_getChrId(gr) (gr)->chrId

char *GenomicRange_setChrName(GenomicRange *gr, char *chrName);

void GenomicRange_free(GenomicRange *range);


#endif
