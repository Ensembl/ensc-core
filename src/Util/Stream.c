#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>

/* for variable argument lists on HP*/
#ifdef __hpux
 #include <varargs.h>
/* for variable argument lists on SGs */
#else
 #include <stdarg.h>
#endif

#define NOEXTERN
#include "Stream.h"
#undef NOEXTERN

#include "Error.h"
#include "StrUtil.h"

#include "FileUtil.h"

STREAM *(InStreams[MAXINSTREAM]);
STREAM  StreamArray[MAXINSTREAM];
int     StreamDepth = 0;
STREAM *OutStream;
STREAM *HistStream;
STREAM *ErrStream;
STREAM *DBGStream;

#ifdef __hpux
/******************************************************************************/
/* Routine    :                                                               */
/*             Stream_fprintf()                                               */
/* Role       :                                                               */
/*             An fprintf using STREAMS                                       */
/* Arguments  :                                                               */
/*             va_alist                                                       */
/* Notes      :                                                               */
/*             This version is for hp which uses varargs for vfprintf()       */
/* History    :                                                               */
/*             09/06/93 SMJS  Initial Implementation for HP                   */
/******************************************************************************/
int Stream_fprintf(va_alist) va_dcl {
  va_list args;
  char *format;
  STREAM *Out;
  int OKFlag = 1;
  char Str[EXTREMELEN];

  va_start(args);
  Out = (STREAM *)va_arg(args,void *);
  format = va_arg(args,char *);

  if (vsnprintf(Str, EXTREMELEN, format, args)<0) {
    Error_write(EVSPRINTF,"Stream_fprintf for hp",ERR_SEVERE,NULL);
    OKFlag=0;
  }

  if (OKFlag) {
    if (!Stream_outputString(Out,Str)) {
      fprintf(stdout,"Trace: Stream_fprintf for HP\n",Out->Fname);
      OKFlag=0;
    }
  }
  va_end(args);
  return OKFlag;
}

#else

/******************************************************************************/
/* Routine    :                                                               */
/*             Stream_fprintf()                                               */
/* Role       :                                                               */
/*             An fprintf using STREAMS                                       */
/* Arguments  :                                                               */
/*             Out   - the STREAM to output to                                */
/*             format- the format string                                      */
/* Notes      :                                                               */
/*             This version is for sgi which uses stdarg for vfprintf()       */
/* History    :                                                               */
/*             09/06/93 SMJS  Initial Implementation for SG                   */
/******************************************************************************/
int Stream_fprintf(STREAM *Out, char *format, ...) {
  va_list args;
  int OKFlag = 1;
  char Str[EXTREMELEN];

  if (Out->File==NULL) {
    printf("Caught Out->File error\n");
    return 0;
  }
  va_start(args,format);

  if (vsnprintf(Str, EXTREMELEN, format, args)<0) {
    Error_write(EVSPRINTF,"Stream_fprintf for SG",ERR_SEVERE,NULL);
    OKFlag=0;
  }

  if (OKFlag) {
    if (!Stream_outputString(Out,Str)) {
      fprintf(stdout,"Trace: Stream_fprintf for SG FName = %s\n",Out->Fname);
      OKFlag=0;
    }
  }

  va_end(args);
  return OKFlag;
}
#endif

/******************************************************************************/
/* Routine    :                                                               */
/*             Stream_flush()                                                 */
/* Role       :                                                               */
/*             To flush a STREAM                                              */
/* Arguments  :                                                               */
/*             Stream - the STREAM                                            */
/* Returns    :                                                               */
/*             None                                                           */
/* History    :                                                               */
/*             09/06/93 SMJS  Initial Implementation                          */
/******************************************************************************/
void Stream_flush(STREAM *Stream) {
  if (Stream->Type==UNIXSTREAM) {
    fflush(Stream->File);
  } else {
    Error_write(EUNKSTYPE,"Stream_flush",ERR_SEVERE,Stream->Fname);
  }
  return;
}

/******************************************************************************/
/* Routine    :                                                               */
/*             Stream_outputString()                                          */
/* Role       :                                                               */
/*             To output a string to a STREAM.                                */
/* Arguments  :                                                               */
/*             Out - the STREAM to output to                                  */
/*             Str - the string to output                                     */
/* Returns    :                                                               */
/*             1 - success                                                    */
/*             0 - failure                                                    */
/* History    :                                                               */
/*             09/06/93 SMJS  Initial Implementation                          */
/******************************************************************************/
int Stream_outputString(STREAM *Out, char *Str) {

  if (Out->Type==UNIXSTREAM) {
    if (Out->File==NULL) {

      fprintf(stderr,"Out->File = NULL\n");
      fprintf(stderr,"Stream_outputString() failed to output string\n");  
      fprintf(stderr,"Str = |%s| ErrNo = %d\n",Str,errno);

    } else if (fprintf(Out->File,"%s",Str)<0) {

      fprintf(stderr,"Stream_outputString() failed to output string\n");
      fprintf(stderr,"Str = |%s| ErrNo = %d\n",Str,errno);
      errno=0;
      clearerr(Out->File);
      return 0;

    } 

    if (Out->TeeFile!=NULL) {
      if (fprintf(Out->TeeFile,"%s",Str)<0) {

        fprintf(stderr,"Stream_outputString() failed to output string to teefile\n");
        fprintf(stderr,"Str = |%s|\n",Str);
        return 0;

      }
    }
  } else {
    fprintf(stderr,"Unknown STREAM Type. Out->Type = %d\n",Out->Type);
  }
  return 1;
}

/******************************************************************************/
/* Routine    :                                                               */
/*             Stream_setDefaults()                                           */
/* Role       :                                                               */
/*             To set the default streams                                     */
/* Arguments  :                                                               */
/*             HistFlag - whether or not to open the Hist.log file            */
/* Returns    :                                                               */
/*             1 - success                                                    */
/*             0 - failure                                                    */
/* History    :                                                               */
/*             09/06/93 SMJS  Initial Implementation                          */
/*             12/05/94 SMJS  Added FilNo to Stream structure for graphics    */
/*             13/02/97 SMJS  Added argument so that we can switch off history*/
/*                            file opening.                                   */
/*             04/05/97 SMJS  Changed to use array of InStreams               */
/******************************************************************************/
int Stream_setDefaults(int HistFlag) {
  char *DBGName;
  int i;

  for (i=0;i<MAXINSTREAM;i++) {
    InStreams[i]=&(StreamArray[i]);
  }

  OutStream = (STREAM *)calloc(1,sizeof(STREAM));
  ErrStream = (STREAM *)calloc(1,sizeof(STREAM));
  HistStream = (STREAM *)calloc(1,sizeof(STREAM));
  DBGStream = (STREAM *)calloc(1,sizeof(STREAM));

  if (OutStream==NULL || ErrStream==NULL || 
      HistStream==NULL || DBGStream==NULL) {
/* Note ErrStream has not been setup */
    fprintf(stderr,"Routine: Stream_setDefaults()\n");
    fprintf(stderr,"Failed to malloc one or more of :-\n");
    fprintf(stderr,"InStream,OutStream,ErrStream,HistStream,DBGStream\n");
    return 0;
  }

  OutStream->File = NULL;
  ErrStream->File = NULL;
  HistStream->File = NULL;
  DBGStream->File = NULL;

/* Setup the Streams */
  InStreams[0]->Type = UNIXSTREAM;
  InStreams[0]->Mode = INTERACTIVE;
  strcpy(InStreams[0]->Fname, "stdin");
  InStreams[0]->File = stdin;
  InStreams[0]->TeeFile = NULL;
  InStreams[0]->FilNo=fileno(stdin);

/* Output to stdout */
  OutStream->Type = UNIXSTREAM;
  OutStream->Mode = INTERACTIVE;
  strcpy(OutStream->Fname, "stdout");
  OutStream->File = stdout;
  OutStream->TeeFile = NULL;
  OutStream->FilNo=fileno(stdout);

/* Debug to stdout */
  DBGStream->Type = UNIXSTREAM;
  DBGStream->Mode = INTERACTIVE;
  strcpy(DBGStream->Fname, "stdout");
  DBGStream->File = stdout;
  DBGStream->TeeFile = NULL;
  DBGStream->FilNo=fileno(stdout);

/* Error to stdout */
  ErrStream->Type = UNIXSTREAM;
  ErrStream->Mode = INTERACTIVE;
  strcpy(ErrStream->Fname, "stdout");
  ErrStream->File = stdout;
  ErrStream->TeeFile = NULL;
  ErrStream->FilNo=fileno(stdout);

/* History to temporary history file */
  HistStream->Type = UNIXSTREAM;
  HistStream->Mode = INTERACTIVE;
  if (HistFlag) {
    strcpy(HistStream->Fname, "Hist.log");
    if ((HistStream->File = FileUtil_open("Hist.log","w","Stream_setDefaults"))==NULL) {
      Error_trace("Stream_setDefaults",HistStream->Fname);
      return 0;
    }
    HistStream->TeeFile = NULL;
    HistStream->FilNo=fileno(HistStream->File);
  }

  return 1;
}

/******************************************************************************/
/* Routine    :                                                               */
/*             Stream_setFile()                                               */
/* Role       :                                                               */
/*             To set the file for a specified output stream.                 */
/*             Two special files are stdout and stderr which are used for the */
/*             standard error and output streams. stdin is not allowed for an */
/*             output stream.                                                 */
/* Arguments  :                                                               */
/*             StreamType - the stream type string                            */
/*             FileName   - the file to set the stream to                     */
/* Returns    :                                                               */
/*             1 - success                                                    */
/*             0 - failure                                                    */
/* History    :                                                               */
/*             26/08/93 SMJS  Initial Implementation                          */
/*             12/05/94 SMJS  Added FilNo to Stream type for graphics         */
/******************************************************************************/
int Stream_setFile(char *StreamType, char *FileName) {
  STREAM *Stream;

  if (!strcmp(StreamType,"OutStream")) {
    Stream=OutStream;
  } else if (!strcmp(StreamType,"DBGStream")) {
    Stream=DBGStream;
  } else if (!strcmp(StreamType,"ErrStream")) {
    Stream=ErrStream;
  } else {
    Error_write(EUNKNOWNSTREAM,"Stream_setFile",ERR_SEVERE,StreamType);
    return 0;
  }

  if (!strcmp(FileName,"stdin")) {
    Error_write(EINVALIDSTREAM,"Stream_setFile",ERR_SEVERE,"stdin not allowed for output");
    return 0;
  }
  if (strcmp(Stream->Fname,"stdout") && 
      strcmp(Stream->Fname,"stderr")) {
    FileUtil_close(Stream->File,Stream->Fname);
  }
  if (!strcmp(FileName,"stdout")) {
    Stream->File = stdout;
  } else if (!strcmp(FileName,"stderr")) {
    Stream->File = stderr;
  } else {
    if ((Stream->File = FileUtil_open(FileName,"w","Stream_setFile"))==NULL) {
      Error_trace("Stream_setFile",FileName);
      return 0;
    }
  }
  Stream->Type = UNIXSTREAM;
  Stream->Mode = INTERACTIVE;
  strcpy(Stream->Fname, FileName);
  Stream->FilNo=fileno(Stream->File);
  return 1;
}

/******************************************************************************/
/* Routine    :                                                               */
/*             Stream_setTeeFile()                                            */
/* Role       :                                                               */
/*             To set the tee file for a specified stream.                    */
/*             Two special files are stdout and stderr which are used for the */
/*             standard error and output streams. stdin is not allowed for an */
/*             output stream.                                                 */
/* Arguments  :                                                               */
/*             StreamType - the stream type string                            */
/*             FileName   - the file to set the stream to                     */
/* Returns    :                                                               */
/*             1 - success                                                    */
/*             0 - failure                                                    */
/* History    :                                                               */
/*             26/08/93 SMJS  Initial Implementation                          */
/******************************************************************************/
int Stream_setTeeFile(char *StreamType, char *FileName) {
  STREAM *Stream;

  if (!strcmp(StreamType,"OutStream")) {
    Stream=OutStream;
  } else if (!strcmp(StreamType,"DBGStream")) {
    Stream=DBGStream;
  } else if (!strcmp(StreamType,"ErrStream")) {
    Stream=ErrStream;
  } else {
    Error_write(EUNKNOWNSTREAM,"Stream_setTeeFile",ERR_SEVERE,StreamType);
    return 0;
  }

  if (Stream->TeeFile!=NULL) {
    fclose(Stream->TeeFile);
  }
  Stream->TeeFile=NULL;
  if (!strcmp(FileName,"stdin") ||
      !strcmp(FileName,"stdout") ||
      !strcmp(FileName,"stderr")) {
    Error_write(EINVALIDSTREAM,"Stream_setTeeFile",ERR_SEVERE,"standard streams not allowed");
    return 0;
  }
  if ((Stream->TeeFile = FileUtil_open(FileName,"w","Stream_setTeeFile"))==NULL) {
    Error_trace("Stream_setTeeFile",FileName);
    return 0;
  }

  return 1;
}
