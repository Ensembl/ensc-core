/*
 * Copyright [1999-2017] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/* for variable argument lists on HP*/
#ifdef __hpux
#include <varargs.h>
/* for variable argument lists on SGs */
#else
#include <stdarg.h>
#endif


#include "Message.h"

#define NOEXTERN
#include "Error.h"
#undef NOEXTERN

/* for ErrStream */
#include "Stream.h"
/* for EXTREMELEN */
#include "StrUtil.h"

int ExitLevel = ERR_FATAL;
int LastErrNo = 0;
int LastELev = 0;

#ifdef __hpux
/******************************************************************************/
/* Routine    :                                                               */
/*             Error_write()                                                  */
/* Role       :                                                               */
/*             To write error messages                                        */
/* Arguments  :                                                               */
/*             va_alist - the variable argument list for this routine. The    */
/*                        first three arguments are :-                        */
/*                        the error number                                    */
/*                        the routine in which the error occurred             */
/*                        a format string for the extra text                  */
/*                        The other arguments are used for the extra string   */
/* Returns    :                                                               */
/*             None                                                           */
/* History    :                                                               */
/*             ??/??/93 SMJS  Initial Implementation                          */
/*             18/10/93 SMJS  Added variable arguments so that the extra      */
/*                            string can say something useful!                */
/*             25/03/98 SMJS  Added error level                               */
/******************************************************************************/
void Error_write(va_alist) va_dcl {
  va_list args;
  char    Str[EXTREMELEN];
  char   *routine;
  char   *extra;
  int     errlev;
  int     errorno;

  va_start(args);
  errorno = va_arg(args,int);
  routine = va_arg(args,char *);
  errlev = va_arg(args,int);
  extra = va_arg(args,char *);

  LastErrNo = errorno;
  LastELev  = errlev;

  Stream_fprintf(ErrStream,"Routine %s\n",routine);

  if (extra != NULL) {
    if (vsnprintf(Str, EXTREMELEN, extra, args)<0) {
      Error_write(EVSPRINTF,"Stream_fprintf for HP",ERR_SEVERE,NULL);
      return;
    }
    Stream_fprintf(ErrStream,"%s :%s\n",message[errorno],Str);

  } else {
    Stream_fprintf(ErrStream,"%s\n",message[errorno]);
  }
  Stream_flush(ErrStream);
  va_end(args);

  if (errlev>=ExitLevel) {
    Stream_fprintf(ErrStream, 
            "Exiting because errlev (%d) is equal to or greater than ExitLevel (%d)\n",
            errlev,ExitLevel);
    exit(1);
  }
}

#else

/******************************************************************************/
/* Routine    :                                                               */
/*             Error_write() (for SG, Linux, Mac, DOS etc)                    */
/* Role       :                                                               */
/*             To write error messages                                        */
/* Arguments  :                                                               */
/*             errorno - the error number                                     */
/*             routine - the routine Error_write was called from              */
/*             extra   - some extra text to add after the standard error text */
/*             ...     - the arguments for the extra string (varargs)         */
/* Returns    :                                                               */
/*             None                                                           */
/* History    :                                                               */
/*             ??/??/93 SMJS  Initial Implementation                          */
/*             18/10/93 SMJS  Added variable arguments so that the extra      */
/*                            string can say something useful!                */
/*             25/03/98 SMJS  Added error level                               */
/******************************************************************************/
void Error_write(int errorno,char *routine,int errlev,char *extra,...) {
  va_list args;
  char    Str[EXTREMELEN];

  va_start(args,extra);

  LastErrNo = errorno;
  LastELev  = errlev;

  Stream_fprintf(ErrStream,"Routine %s\n",routine);

  if (extra != NULL) {
    if (vsprintf(Str, extra, args)<0) {
      Error_write(EVSPRINTF,"Stream_fprintf for SG",ERR_SEVERE,NULL);
      return;
    }
    Stream_fprintf(ErrStream,"%s :%s\n",message[errorno],Str);
  } else {
    Stream_fprintf(ErrStream,"%s\n",message[errorno]);
  }
  Stream_flush(ErrStream);
  va_end(args);

  if (errlev>=ExitLevel) {
    Stream_fprintf(ErrStream, 
            "Exiting because errlev (%d) is equal to or greater than ExitLevel (%d)\n",
            errlev,ExitLevel);
    exit(1);
  }
}
#endif

/******************************************************************************/
/* Routine    :                                                               */
/*             Error_trace()                                                  */
/* Role       :                                                               */
/*             To print a message indicating that an error was returned from  */
/*             a routine, indicating the current routine (for tracing)        */
/* Arguments  :                                                               */
/*             Routine - the name of the current routine                      */
/*             Info    - an information string (can be NULL)                  */
/* Returns    :                                                               */
/*             None                                                           */
/******************************************************************************/
void Error_trace(char *Routine, char *Info) {
  Stream_fprintf(ErrStream,"Trace: Error returned to %s().\n",Routine);
  if (Info!=NULL) {
    Stream_fprintf(ErrStream,"Message: %s\n",Info);
  }
  Stream_flush(ErrStream);
}
