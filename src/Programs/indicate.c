/****************************************************************\
*                                                                *
*  indicate : Index a file using OBDA indexers                   *
*                                                                *
*  Steve Searle.        mailto:searle@sanger.ac.uk               *
*  Copyright (C) 2004.  All Rights Reserved.                     *
*                                                                *
*  This source code is distributed under the terms of the        *
*  GNU Lesser General Public License. See the file COPYING       *
*  or http://www.fsf.org/copyleft/lesser.html for details        *
*                                                                *
*  If you use this code, please keep this notice intact.         *
*                                                                *
\****************************************************************/

#include <stdio.h>
#include <stdlib.h> /* For atoi() */
#include <sys/types.h>
#include <sys/stat.h>
#include <libgen.h>
#include <string.h>

#include "indicate.h"
#include "BioIndex.h"
#include "parsers.h"
#include "Error.h"

void Indicate_usage();

int main(int argc, char **argv){
  char *arg;
  char *file_prefix = NULL;
  char *index_dir = NULL;
  char *data_dir = NULL;
  char *primary_parser = NULL;
  char *secondary_parser = NULL;
  char *secondary_multi_parser = NULL;
  BioIndex_Parser_Func parser;
  BioIndex_Parser_Func sec_parser;
  BioIndex_MultiParser_Func sec_multi_parser;
  char *line_prefix = NULL;
  char *secondary_line_prefix = NULL;
  char *record_prefix = NULL;
  BioIndex_Index_Definition *primary_def;
  BioIndex_Index_Definition *secondary_def;
  Vector *secondary_defs = Vector_new();
  BioIndex *bi;
  struct stat dirStat;
  char *gtStr = ">";
  char *emptyStr = "";
  int argNum = 1;

  initEnsC();

  primary_parser = secondary_parser = secondary_multi_parser = emptyStr;
  line_prefix = secondary_line_prefix = record_prefix = gtStr;

  while (argNum < argc) {
    char *arg = argv[argNum];
    char *val;
    
    if (argNum == argc-1) {
      Indicate_usage();
    }

    val = argv[++argNum];
//    printf("%s %s\n",arg,val);

    if (!strcmp(arg, "-d") || !strcmp(arg,"--data_dir")) {
      StrUtil_copyString(&data_dir,val,0);
    } else if (!strcmp(arg, "-f") || !strcmp(arg,"--file_prefix")) {
      StrUtil_copyString(&file_prefix,val,0);
    } else if (!strcmp(arg, "-i") || !strcmp(arg,"--index")) {
      StrUtil_copyString(&index_dir,val,0);
    } else if (!strcmp(arg, "-p") || !strcmp(arg,"--parser")) {
      StrUtil_copyString(&primary_parser,val,0);
    } else if (!strcmp(arg, "-s") || !strcmp(arg,"--secondary_parser")) {
      StrUtil_copyString(&secondary_parser,val,0);
    } else if (!strcmp(arg, "-M") || !strcmp(arg,"--secondary_multi_parser")) {
      StrUtil_copyString(&secondary_multi_parser,val,0);
    } else if (!strcmp(arg, "-m") || !strcmp(arg,"--secondary_line_prefix")) {
      StrUtil_copyString(&secondary_line_prefix,val,0);
    } else if (!strcmp(arg, "-l") || !strcmp(arg,"--line_prefix")) {
      StrUtil_copyString(&line_prefix,val,0);
    } else if (!strcmp(arg, "-r") || !strcmp(arg,"--record_prefix")) {
      StrUtil_copyString(&record_prefix,val,0);
    }

    argNum++;
  }

  printf("Program for generating OBDA indices\n"
         "Steve M.J. Searle.  searle@sanger.ac.uk  Last update Feb 2004.\n");

  if (!file_prefix || !index_dir) {
    Indicate_usage();
  }

  printf("Using File:  [%s]\n"
         "      Index: [%s]\n", file_prefix, index_dir);

  if (!stat(index_dir,&dirStat)) {
    if (!(dirStat.st_mode & S_IFDIR)) {
      Error_write(EOBDA,"main",ERR_SEVERE,"Output dir %s is not a directory",index_dir);
    }
  } else {
    printf("Creating Index dir %s", index_dir);
    mkdir(index_dir,0755);
  }

  printf("Primary parser = %s\n",primary_parser);
  parser = Parser_lookup(primary_parser);
  printf("Got function = %d\n",parser);
  primary_def = BioIndex_Index_Definition_create(record_prefix,
                                                 line_prefix,
                                                 parser,
                                                 NULL,
                                                 "ACC");
  if (strlen(secondary_multi_parser)) {
    printf("Secondary multi parser = %s\n",secondary_multi_parser);
    sec_multi_parser = Parser_lookup(secondary_multi_parser);
    printf("Got function = %d\n",secondary_parser);
    secondary_def = BioIndex_Index_Definition_create(record_prefix,
                                                     secondary_line_prefix,
                                                     NULL,
                                                     sec_multi_parser,
                                                     "ID");
    Vector_addElement(secondary_defs,secondary_def);
  }

  if (strlen(secondary_parser)) {
    printf("Secondary parser = %s\n",secondary_parser);
    sec_parser = Parser_lookup(secondary_parser);
    printf("Got function = %d\n",sec_parser);
    secondary_def = BioIndex_Index_Definition_create(record_prefix,
                                                     secondary_line_prefix,
                                                     sec_parser,
                                                     NULL,
                                                     "ID");
    Vector_addElement(secondary_defs,secondary_def);
  }
 
  bi = BioIndex_generate_flat(index_dir,data_dir,file_prefix,data_dir,primary_def,
                              secondary_defs); 

  printf("Completed Indexing");
  return 0;
}

void Indicate_usage() {
  printf("indicate \n"
         "-d --datadir Path to files to be indexed\n"
         "-f --file_prefix File prefix of files to be indexed\n"
         "-i --index Output index dir name\n"
         "-p --parser Parser function used to parse the primary key line\n"
         "-s --secondary_parser Parser function used to parse the secondary key line\n"
         "-M --secondary_multi_parser Parser function used to parse the secondary key line\n"
         "-m --secondary_line_prefix The string used to identify secondary ID lines\n"
         "-l --line_prefix The string used to identify primary ID lines\n"
         "-r --record_prefix Line prefix at start of record\n");
  exit(1);
}
