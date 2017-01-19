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

#ifndef __STREAM_H__
#define __STREAM_H__

#include <sys/param.h>

/* Modes for input */
#define INTERACTIVE 1               /* Interactive input mode */
#define COMFILE     2               /* Comfile input mode */
#define SOCKET      3               /* Socket */

#define MAXINSTREAM 100

/* Stream Types (only one)  */
#define UNIXSTREAM  1               /* Standard unix stream */

typedef struct Stream STREAM;

/* Structure to define a stream for input or output */
struct Stream {
  int   Mode;             /* Mode for file (COMFILE or INTERACTIVE) */
  int   Type;             /* the type of the stream */
  FILE *File;             /* File pointer to stream */
  char  Fname[MAXPATHLEN]; /* the filename */
  FILE *TeeFile;          /* A file to also send interactive data */
  int   FilNo;            /* File descriptor number for stream */
};

#ifndef NONEXTERN
  extern STREAM *(InStreams[MAXINSTREAM]);
  extern int     StreamDepth;
  extern STREAM *OutStream;
  extern STREAM *ErrStream;
  extern STREAM *HistStream;
  extern STREAM *DBGStream;
#endif /* NOEXTERN */

/* prototypes for stream */
#ifdef __hpux
 int       Stream_fprintf();
#else
 int       Stream_fprintf(STREAM *Out, char *format, ...);
#endif
void      Stream_flush(STREAM *Stream);
int       Stream_outputString(STREAM *Out, char *Str);
int       Stream_setDefaults(int HistFlag);
int       Stream_setFile(char *StreamType, char *FileName);
int       Stream_setTeeFile(char *StreamType, char *FileName);

#endif /* __STREAM_H__ */
