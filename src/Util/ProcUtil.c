#include "ProcUtil.h"

#include <stdio.h>
#include <string.h>
#ifdef linux
#include <malloc.h>
#endif

#ifndef dos
#include <unistd.h>
#include <sys/times.h>
#endif

#include "FileUtil.h"
/* for error handling */
#include "Error.h"
/* for OutStream */
#include "Stream.h"
/******************************************************************************/
/* Routine    :                                                               */
/*             TimeInfo()                                                     */
/* Role       :                                                               */
/*             timing information about a process                             */
/* Arguments  :                                                               */
/*             routine - the routine calling                                  */
/* Returns    :                                                               */
/*             None                                                           */
/******************************************************************************/
int ProcUtil_timeInfo(char *routine) {
#ifndef dos
  struct tms t;
  int CLOK_TCK;

  CLOK_TCK = (int)sysconf(_SC_CLK_TCK);

  times(&t);

  FPrintF(OutStream,"----------------------------------\n");
  if (strlen(routine)) {
    FPrintF(OutStream,"Timing information %s\n",routine);
  }
  FPrintF(OutStream,"User time      : %f s\n",(float)t.tms_utime/(float)CLOK_TCK);
  FPrintF(OutStream,"System time    : %f s\n",(float)t.tms_stime/(float)CLOK_TCK);
  FPrintF(OutStream,"User children  : %f s\n",(float)t.tms_cutime/(float)CLOK_TCK);
  FPrintF(OutStream,"System children: %f s\n",(float)t.tms_cstime/(float)CLOK_TCK);
  FPrintF(OutStream,"----------------------------------\n");
#endif

  return 1;
}

/******************************************************************************/
/* Routine   :                                                                */
/*            mallInfo()                                                      */
/* Role      :                                                                */
/*            To output the memory allocation information from                */
/*            the mallinfo structure                                          */
/* Arguments :                                                                */
/*            None                                                            */
/* Returns   :                                                                */
/*            None                                                            */
/* History   :                                                                */
/*            22/6/93 SMJS  Initial Implementation                            */
/******************************************************************************/
void ProcUtil_mallInfo() {
#if !defined(os2) && !defined(__ppc__)
#ifndef linux
  struct mallinfo mall;

  mall=mallinfo();

  FPrintF(OutStream,"---------------------------------------------------\n\n");
  FPrintF(OutStream,"Total space in arena           : %d\n", mall.arena);
  FPrintF(OutStream,"Number of ordinary blocks      : %d\n", mall.ordblks);
  FPrintF(OutStream,"Number of small blocks         : %d\n", mall.smblks);
  FPrintF(OutStream,"Number of holding blocks       : %d\n", mall.hblks);
  FPrintF(OutStream,"Space in holding block headers : %d\n", mall.hblkhd);
  FPrintF(OutStream,"Space in small blocks in use   : %d\n", mall.usmblks);
  FPrintF(OutStream,"Space in free small blocks     : %d\n", mall.fsmblks);
  FPrintF(OutStream,"Space in ordinary blocks in use: %d\n", mall.uordblks);
  FPrintF(OutStream,"Space in free ordinary blocks  : %d\n", mall.fordblks);
  FPrintF(OutStream,"Cost of enabling keep option   : %d\n", mall.keepcost);
  FPrintF(OutStream,"---------------------------------------------------\n\n");
#endif

#ifdef linux_old
#ifdef linux_orig
  struct mstats stats;
  stats=mstats();
  FPrintF(OutStream,"---------------------------------------------------\n\n");
  FPrintF(OutStream,"Total space in heap            : %d\n",stats.bytes_total);
  FPrintF(OutStream,"Chunks allocated               : %d\n",stats.chunks_used);
  FPrintF(OutStream,"Bytes used                     : %d\n",stats.bytes_used);
  FPrintF(OutStream,"Chunks free                    : %d\n",stats.chunks_free);
  FPrintF(OutStream,"Bytes free                     : %d\n",stats.bytes_free);
  FPrintF(OutStream,"---------------------------------------------------\n\n");
#endif
#endif
#endif
}

