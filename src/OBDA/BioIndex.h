/****************************************************************\
*                                                                *
*  Library for manipulation of Bio format indexes                *
*                                                                *
*  Steve Searle         mailto:searle@sanger.ac.uk               *
*  Copyright (C) 2002.  All Rights Reserved.                     *
*                                                                *
*  This source code is distributed under the terms of the        *
*  GNU Lesser General Public License. See the file COPYING       *
*  or http://www.fsf.org/copyleft/lesser.html for details        *
*                                                                *
*  If you use this code, please keep this notice intact.         *
*                                                                *
\****************************************************************/

#ifndef INCLUDED_BIOINDEX_H
#define INCLUDED_BIOINDEX_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "Vector.h"

#include <stdio.h>

#include "Sequence.h"


typedef enum {
  BioIndexType_FLAT    = (1<<0),
  BioIndexType_BDB     = (1<<1)
} BioIndexType;  

typedef struct {
  int width;
  char type;
} BioIndex_FieldProperties;

typedef struct {
  void *handle;
  char *file_name;
  int   rec_length;
  int   num_record;
  int   num_field;
  BioIndex_FieldProperties  *field_properties;
  long  header_len;
  char *namespace;
} BioIndex_IndexFile; 

typedef struct {
  BioIndex_IndexFile *primary_index;
  BioIndex_IndexFile **secondary_indices;
  int       num_secondary;
  char    **secondary_namespaces;
  char     *path;
  Vector *fileids;
  char     *primary_namespace;
} BioIndex;

typedef struct {
  FILE   *fp;
  int    id;
  char  *name;
  off_t   length;
} BioIndex_FileID;

typedef struct{
  long   start;
  int    length;
  char *entry;
  short  fileIndex;
} BioIndex_Location;

typedef char * (*BioIndex_Parser_Func)(char *line);
typedef char * (*BioIndex_MultiParser_Func)(char *line,char **retPos);
typedef struct {
  char *record_start;
  char *line_prefix;
  char *type;
  int   prefix_len;
  int   start_len;
  BioIndex_Parser_Func parser;
  BioIndex_MultiParser_Func multiParser;
} BioIndex_Index_Definition; 

BioIndex  *BioIndex_open(char *path);
    void   BioIndex_close(BioIndex *bi);
    void   BioIndex_rewind(BioIndex *bi);

Sequence **BioIndex_fetch_all(BioIndex *bi, char *key);
Sequence *BioIndex_fetch(BioIndex *bi, char *key);
Sequence * BioIndex_get_fasta_sequence(BioIndex *bi, BioIndex_Location *loc,
                                       char *id);

char * BioIndex_get_entry(BioIndex *bi, BioIndex_Location *loc);
BioIndex_Index_Definition *BioIndex_Index_Definition_create(char *record_start,
                                              char *line_prefix,
                                              BioIndex_Parser_Func parser,
                                              BioIndex_MultiParser_Func multi_parser,
                                              char *type);
Vector * BioIndex_get_by_secondary_key(BioIndex *bi, char *key,
                                          BioIndex_IndexFile *index);
BioIndex *BioIndex_open_flat(char *path);

BioIndex *BioIndex_generate_flat(char *path, char *seq_path,
                                 char *select, char *format,
                                 BioIndex_Index_Definition *primary_def,
                                 Vector *secondary_defs);
#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_BIOINDEX_H */
