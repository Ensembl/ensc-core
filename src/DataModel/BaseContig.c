#define __BASECONTIG_MAIN__
#include "BaseContig.h"
#undef __BASECONTIG_MAIN__

void BaseContig_freePtrs(BaseContig *bc) {
  Sequence_freePtrs((Sequence *)bc);
}
