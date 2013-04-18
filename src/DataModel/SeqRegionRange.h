#ifndef __SEQREGIONRANGE_H__
#define __SEQREGIONRANGE_H__

#include "DataModelTypes.h"
#include "EnsC.h"

struct SeqRegionRangeStruct {
  char * srName;
  int    srStart;
  int    srEnd;
  IDType srId;
};

SeqRegionRange *SeqRegionRange_new();

#define SeqRegionRange_setSeqRegionStart(srr,s) (srr)->srStart = (s)
#define SeqRegionRange_getSeqRegionStart(srr) (srr)->srStart

#define SeqRegionRange_setSeqRegionEnd(srr,e) (srr)->srEnd = (e)
#define SeqRegionRange_getSeqRegionEnd(srr) (srr)->srEnd

char *setSeqRegionName(SeqRegionRange *srRange, char *srName);
#define SeqRegionRange_getSeqRegionName(srr) (srr)->srName

#define SeqRegionRange_setSeqRegionId(srr,i) (srr)->srId = (i)
#define SeqRegionRange_getSeqRegionId(srr) (srr)->srId

char *SeqRegionRange_setSeqRegionName(SeqRegionRange *srRange, char *srName);

void SeqRegionRange_free(SeqRegionRange *range);

void SeqRegionRange_expand(SeqRegionRange *range, int pad);


#endif
