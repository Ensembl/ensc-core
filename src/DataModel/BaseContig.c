#define __BASECONTIG_MAIN__
#include "BaseContig.h"
#undef __BASECONTIG_MAIN__
#include "Slice.h"
#include "RawContig.h"

/* Needs to do polymorphic call as Slice name is auto generated */
#ifdef OLD
char *BaseContig_getName(BaseContig *bc) {
  switch (BaseContig_getObjectType(bc)) {
    case CLASS_RAWCONTIG:
      return RawContig_getName((RawContig *)bc);
      break;
    case CLASS_SLICE:
      return Slice_getName((Slice *)bc);
      break;
    default:
      fprintf(stderr,"ERROR: Unknown contig type in getName\n");
      exit(1);
  }
}
#endif

