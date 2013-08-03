#ifndef __PROCUTIL_H__
#define __PROCUTIL_H__

#ifdef linux
#define UNW_LOCAL_ONLY
//#include "libunwind.h"
#endif

int ProcUtil_timeInfo(char *routine);
void ProcUtil_mallInfo(void);

#if 0
void ProcUtil_showBacktrace(char *prog);
#ifdef linux
int ProcUtil_getFileAndLine(char *prog, unw_word_t addr, char *file, size_t flen, int *line);
#endif
#endif

#endif
