#include "ProcUtil.h"

#include <stdio.h>
#include <string.h>
#ifdef linux
#include <malloc.h>
#endif

#ifdef linux
#define UNW_LOCAL_ONLY
#include "libunwind.h"
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
/*             ProcUtil_timeInfo()                                            */
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

  Stream_fprintf(OutStream,"----------------------------------\n");
  if (strlen(routine)) {
    Stream_fprintf(OutStream,"Timing information %s\n",routine);
  }
  Stream_fprintf(OutStream,"User time      : %f s\n",(float)t.tms_utime/(float)CLOK_TCK);
  Stream_fprintf(OutStream,"System time    : %f s\n",(float)t.tms_stime/(float)CLOK_TCK);
  Stream_fprintf(OutStream,"User children  : %f s\n",(float)t.tms_cutime/(float)CLOK_TCK);
  Stream_fprintf(OutStream,"System children: %f s\n",(float)t.tms_cstime/(float)CLOK_TCK);
  Stream_fprintf(OutStream,"----------------------------------\n");
#endif

  return 1;
}

/******************************************************************************/
/* Routine   :                                                                */
/*            ProcUtil_mallInfo()                                             */
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
#if !defined(os2) && !defined(__ppc__) && !defined(__osf__)
#ifndef linux
  struct mallinfo mall;

  mall=mallinfo();

  Stream_fprintf(OutStream,"---------------------------------------------------\n\n");
  Stream_fprintf(OutStream,"Total space in arena           : %d\n", mall.arena);
  Stream_fprintf(OutStream,"Number of ordinary blocks      : %d\n", mall.ordblks);
  Stream_fprintf(OutStream,"Number of small blocks         : %d\n", mall.smblks);
  Stream_fprintf(OutStream,"Number of holding blocks       : %d\n", mall.hblks);
  Stream_fprintf(OutStream,"Space in holding block headers : %d\n", mall.hblkhd);
  Stream_fprintf(OutStream,"Space in small blocks in use   : %d\n", mall.usmblks);
  Stream_fprintf(OutStream,"Space in free small blocks     : %d\n", mall.fsmblks);
  Stream_fprintf(OutStream,"Space in ordinary blocks in use: %d\n", mall.uordblks);
  Stream_fprintf(OutStream,"Space in free ordinary blocks  : %d\n", mall.fordblks);
  Stream_fprintf(OutStream,"Cost of enabling keep option   : %d\n", mall.keepcost);
  Stream_fprintf(OutStream,"---------------------------------------------------\n\n");
#endif

#ifdef linux_old
#ifdef linux_orig
  struct mstats stats;
  stats=mstats();
  Stream_fprintf(OutStream,"---------------------------------------------------\n\n");
  Stream_fprintf(OutStream,"Total space in heap            : %d\n",stats.bytes_total);
  Stream_fprintf(OutStream,"Chunks allocated               : %d\n",stats.chunks_used);
  Stream_fprintf(OutStream,"Bytes used                     : %d\n",stats.bytes_used);
  Stream_fprintf(OutStream,"Chunks free                    : %d\n",stats.chunks_free);
  Stream_fprintf(OutStream,"Bytes free                     : %d\n",stats.bytes_free);
  Stream_fprintf(OutStream,"---------------------------------------------------\n\n");
#endif
#endif
#endif
}


// Adapted from code from http://blog.bigpixel.ro/2010/09/stack-unwinding-stack-trace-with-gcc/
// Needs libunwind
void ProcUtil_showBacktrace(char *prog) {
#ifdef linux
  char name[2048];
  unw_cursor_t cursor; unw_context_t uc;
  unw_word_t ip, sp, offp;

  unw_getcontext(&uc);
  unw_init_local(&cursor, &uc);

  while (unw_step(&cursor) > 0) {
    char file[2048];
    int line = 0;

    name[0] = '\0';
    unw_get_proc_name(&cursor, name, 2048, &offp);
    unw_get_reg(&cursor, UNW_REG_IP, &ip);
    unw_get_reg(&cursor, UNW_REG_SP, &sp);

    if (prog == NULL) {
      fprintf(stderr, "%s ip = %lx, sp = %lx\n", name, (long)ip, (long)sp);
    } else {
      ProcUtil_getFileAndLine(prog, (long)ip, file, 2048, &line);
      fprintf(stderr, "%s in file %s line %d\n", name, file, line);
    }
  }
#else
  fprintf(stderr, "backtrace not supported on this platform\n");
#endif
}


int ProcUtil_getFileAndLine(char *prog, unw_word_t addr, char *file, size_t flen, int *line) {
#ifdef linux
  static char buf[2048];
  char *p;

  // prepare command to be executed
  // our program need to be passed after the -e parameter
  sprintf (buf, "/usr/bin/addr2line -C -e %s -f -i %lx", prog, addr);
  FILE* f = popen (buf, "r");

  if (f == NULL) {
    perror (buf);
    return 0;
  }

  // get function name
  fgets (buf, 2048, f);

  // get file and line
  fgets (buf, 2048, f);

  if (buf[0] != '?') {
    int l;
    char *p = buf;

    // file name is until ':'
    while (*p != ':') {
      p++;
    }

    *p++ = 0;
    // after file name follows line number
    strcpy (file , buf);
    sscanf (p,"%d", line);
  } else {
    strcpy (file,"unknown");
    *line = 0;
  }

  pclose(f);
#else
  fprintf(stderr, "getFileAndLine not supported on this platform\n");
#endif
}
