#ifndef __COORDPAIR_H__
#define __COORDPAIR_H__

typedef struct CoordPairStruct CoordPair;

struct CoordPairStruct {
  long start;
  long end;
};

CoordPair *CoordPair_new(long start, long end);

#define CoordPair_getEnd(cp) (cp)->end
#define CoordPair_getStart(cp) (cp)->start
#endif
