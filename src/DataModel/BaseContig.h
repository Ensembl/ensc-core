#ifndef __BASECONTIG_H__
#define __BASECONTIG_H__

#include "DataModelTypes.h"

typedef enum ContigTypeEnum {
  CONTIGTYPE_NONE,
  RAWCONTIG,
  SLICE
} ContigType;

#define BASECONTIG_DATA \
  ContigType contigType; \
  int start; \
  int end;

struct BaseContigStruct {
  BASECONTIG_DATA
};

#endif
