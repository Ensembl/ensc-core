#include "BaseContig.h"
#include "Slice.h"
#include "RawContig.h"

/* Needs to do polymorphic call as Slice name is auto generated */
char *BaseContig_getName(BaseContig *bc) {
  switch (BaseContig_getContigType(bc)) {
    case RAWCONTIG:
      return RawContig_getName((RawContig *)bc);
      break;
    case SLICE:
      return Slice_getName((Slice *)bc);
      break;
    default:
      fprintf(stderr,"ERROR: Unknown contig type in getName\n");
      exit(1);
  }
}
