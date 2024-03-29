/*
 * Copyright [1999-2024] EMBL-European Bioinformatics Institute
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *      http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "ProcUtil.h"
#include <stdio.h>
#include "FileUtil.h"
#include "StrUtil.h"
#include <stdlib.h>

#include <stdio.h>
#include <string.h>
#ifdef linux
#include <malloc.h>
#endif

#ifdef linux
#define UNW_LOCAL_ONLY
//#include "libunwind.h"
#endif

#ifdef __APPLE__
#define UNW_LOCAL_ONLY
//#include <libunwind.h>
#include <mach-o/dyld.h>
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

  Stream_fprintf(ErrStream,"----------------------------------\n");
  if (strlen(routine)) {
    Stream_fprintf(ErrStream,"Timing information %s\n",routine);
  }
  Stream_fprintf(ErrStream,"User time      : %f s\n",(float)t.tms_utime/(float)CLOK_TCK);
  Stream_fprintf(ErrStream,"System time    : %f s\n",(float)t.tms_stime/(float)CLOK_TCK);
  Stream_fprintf(ErrStream,"User children  : %f s\n",(float)t.tms_cutime/(float)CLOK_TCK);
  Stream_fprintf(ErrStream,"System children: %f s\n",(float)t.tms_cstime/(float)CLOK_TCK);
  Stream_fprintf(ErrStream,"----------------------------------\n");
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
#ifndef __APPLE__
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

#if 0
int ProcUtil_getFileAndLine(char *prog, unw_word_t addr, char *file, size_t flen, int *line);

// Adapted from code from http://blog.bigpixel.ro/2010/09/stack-unwinding-stack-trace-with-gcc/
// Needs libunwind
void ProcUtil_showBacktrace(char *prog) {
#if defined(linux) || defined(__APPLE__)
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

//#if defined(linux)
    if (prog == NULL) {
      fprintf(stderr, "%s ip = %lx, sp = %lx\n", name, (long)ip, (long)sp);
    } else {
      ProcUtil_getFileAndLine(prog, ip, file, 2048, &line);
      fprintf(stderr, "%s in file %s line %d\n", name, file, line);
    }
//# else
//    fprintf(stderr, "%s ip = %lx, sp = %lx\n", name, (long)ip, (long)sp);
//#endif
  }
#else
  fprintf(stderr, "backtrace not supported on this platform\n");
#endif
}


int ProcUtil_getFileAndLine(char *prog, unw_word_t addr, char *file, size_t flen, int *line) {
  static char buf[2048];
  char *p;

  // prepare command to be executed
  // our program need to be passed after the -e parameter
#if defined(linux)
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
#elif defined(__APPLE__)
/// NOT WORKING YET. Need extra magic to get all the addresses I need
  //sprintf (buf, "/usr/bin/atos -p %s", prog, addr);
  void *mach_header = _dyld_get_image_header(0);
  char exeName[1024];
  FileUtil_stripPath(prog,exeName);
  sprintf (buf, "/usr/bin/atos -p %s -l %p %p", exeName,mach_header,addr);
//  fprintf(stderr,"Command: %s\n",buf);
  FILE* f = popen (buf, "r");

//Examples:
//AttributeAdaptor_storeOnGeneId (in GeneWriteTest) (AttributeAdaptor.c:110)
//GeneAdaptor_store (in GeneWriteTest) + 2060

  if (f == NULL) {
    perror (buf);
    return 0;
  }

/*
  // get function name
  fgets (buf, 2048, f);

  // get file and line
  fgets (buf, 2048, f);
*/

  while (fgets(buf,2048, f)) {
//    fprintf(stderr,"%s\n",buf);
    int nTok;
    char **tokens;
    StrUtil_tokenize(&tokens, &nTok, buf);
    char *colonPos = NULL;
    if ((colonPos = strchr(tokens[nTok-1], ':')) != NULL) {
      int len = colonPos-tokens[nTok-1]-1;
      strncpy(file, &(tokens[nTok-1][1]), len);
      file[len] = '\0';
      sscanf(colonPos+1,"%d",line);
    } else {
      strcpy (file,"unknown");
      *line = 0;
    }
    int i;
    for (i=0;i<nTok;i++) {
      free(tokens[i]);
    }
  }

  pclose(f); 
#else
  fprintf(stderr, "getFileAndLine not supported on this platform\n");
#endif
}
#endif
