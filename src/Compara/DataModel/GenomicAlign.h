#ifndef __GENOMICALIGN_H__
#define __GENOMICALIGN_H__

#include "ComparaDataModelTypes.h"
#include "EnsRoot.h"

struct GenomicAlignStruct {
  ENSROOT_DATA
  Storable st;
  char *cigarString;
  DNAFrag *consensusDNAFrag;
  int consensusStart;
  int consensusEnd;
  DNAFrag *queryDNAFrag;
  int queryStart;
  int queryEnd;
  signed char queryStrand;
  char *alignmentType;
  double percentId;
  double score;
};

#endif
