#ifndef __BASECONTIG_H__
#define __BASECONTIG_H__

typedef struct BaseContigStruct BaseContig;

typedef enum ContigTypeEnum {
  RAWCONTIG,
  SLICE
} ContigType;

#define BASECONTIG_DATA \
  ContigType contigType; \
  int start; \
  int end;

struct BaseContigStruct {
};

#endif
