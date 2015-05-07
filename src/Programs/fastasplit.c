/*
 * Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
//#include <malloc.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>



#define MAXSTRLEN 1024*32
#define HEADER_ALLOC_SIZE 1000

#define _FILE_OFFSET_BITS 64

void printUsage(void);

int main(int argc, char *argv[]) {
  FILE *InP;
  FILE *OutP = NULL;
  char outdir[MAXSTRLEN];
  char line[MAXSTRLEN];
  fpos_t **headerPositions;
  int nHeader = 0;
  fpos_t *currentPos;
  int i;
  long randPos1;
  long randPos2;
  fpos_t *tmpPos;
  fpos_t tmp_pos;   
  long fileSize;
  long currentSize = 0;
  long nChunk;
  char chunkFName[MAXSTRLEN];
  long currentChunk = 0;
  char baseName[MAXSTRLEN];
  char databaseName[MAXSTRLEN];

  int chunkArgNum;
  int dirArgNum;
  int fileArgNum;
  int bySize = 0;
  long nPerFile;
  long entryCount = -1;
  int argnum;
  long remainder; 

  if (argc<4) {
    printUsage();
    exit(1);
  } 

  databaseName[0] = '\0';

  argnum = 1;
  while (argnum < argc-3) {
    char *arg = argv[argnum];
    char *val;
  
    if (!strcmp(arg, "-s") || !strcmp(arg,"--size")) {
      bySize = 1;
    } else if (!strcmp(arg, "-f") || !strcmp(arg,"--file_prefix")) {
      val = argv[++argnum];
      strcpy(databaseName,val);
    } else {
      printUsage();
      exit(1);
    } 

    argnum++;
  }

  if (argc - argnum != 3) {
    printUsage();
    exit(1);
  } 
  
  fileArgNum = argnum++;
  if ((InP = fopen(argv[fileArgNum],"r")) == NULL) {
    fprintf(stderr,"Error: Couldn't open file %s\n",argv[fileArgNum]);
    exit(1);
  }

  if (!strlen(databaseName)) {
    strippath(argv[fileArgNum],baseName); 
    stripext(baseName,databaseName); 
  }

  chunkArgNum = argnum++;
  if (!(nChunk = atol(argv[chunkArgNum]))) {
    fprintf(stderr,"Error: Couldn't parse nChunk into a number %s\n",argv[chunkArgNum]);
    exit(1);
  } 


  dirArgNum = argnum++;
  strcpy(outdir,argv[dirArgNum]);


  if ((headerPositions = (fpos_t **)calloc(HEADER_ALLOC_SIZE,sizeof(fpos_t *)))==NULL) {
    fprintf(stderr,"Error: Allocating headerPositions array\n");
    exit(1);
  }


  fgetpos(InP,&tmp_pos);
  while (fgets(line,MAXSTRLEN,InP)) {
    if (line[0] == '>') {
      fpos_t *pos;
      pos = calloc(1,sizeof(fpos_t));
      memcpy(pos,&tmp_pos,sizeof(fpos_t));
      headerPositions[nHeader++] = pos;
      if (nHeader && (nHeader%HEADER_ALLOC_SIZE) == 0) {
        if ((headerPositions = (fpos_t **)realloc(headerPositions,(nHeader+HEADER_ALLOC_SIZE)*sizeof(fpos_t *)))==NULL) {
          fprintf(stderr,"Error: Allocating headerPositions array\n");
          exit(1);
        }
      }
    }
    fgetpos(InP,&tmp_pos);
  }
  if (bySize) { 
    struct stat infile_stats;
    stat(argv[fileArgNum],&infile_stats);
    fileSize = infile_stats.st_size/nChunk;

  
    if (fileSize < 100) {
      fprintf(stderr,
           "Do you really want to split you're file into pieces of %ld bytes\n",
           fileSize);
      exit(1);
    }
  } else {
    nPerFile = nHeader/nChunk;
    printf("nPerFile = %ld nHeader = %d nChunk = %ld\n",nPerFile,nHeader,nChunk);
    if (nPerFile < 1) {
      fprintf(stderr,
           "Do you really want to split you're file into pieces of with less than 1 sequence per file\n");
      exit(1);
    } else {
      remainder = nHeader - nChunk*nPerFile;
      printf ("Remainder = %ld\n",nHeader - nChunk*nPerFile);
    }
  }


  rewind(InP);

  for (i=0;i<10*nHeader;i++) {
    randPos1 = lrand48()%nHeader;
    randPos2 = lrand48()%nHeader;
    tmpPos = headerPositions[randPos1];
    headerPositions[randPos1] = headerPositions[randPos2];
    headerPositions[randPos2] = tmpPos;
  }
  printf("AFTER RANDOMISATION\n");
  
  currentSize = -1;
  if (remainder) {
    nPerFile++;
  } 
  
  for (i=0;i<nHeader;i++){
    if ((bySize && (currentSize < 0 ||  currentSize  > fileSize)) ||
        (entryCount < 0 || entryCount >= nPerFile)) {
      if (OutP) fclose(OutP);
      if (bySize) {
        currentSize = 0;
      } else {
        entryCount = 0;
      }

      sprintf(chunkFName,"%s/%s_chunk_%7.7ld",
              outdir,databaseName,currentChunk++);

      if (remainder && currentChunk == remainder+1) {
        nPerFile--;
      }

      if ((OutP = fopen(chunkFName,"w")) == NULL) {
        fprintf(stderr,"Error: Couldn't open file %s\n",chunkFName);
        exit(1);
      }
    }
    fsetpos(InP,headerPositions[i]);
    fgets(line,MAXSTRLEN,InP);
    if (bySize) {
      currentSize += strlen(line);
    } else {
      entryCount++;
    }

    fprintf(OutP,"%s",line);
    while (fgets(line,MAXSTRLEN,InP)) {
      if (line[0] == '>') break;
      if (bySize) {
        currentSize += strlen(line);
      }
      fprintf(OutP,"%s",line);
    }  
  }
  if (OutP) fclose(OutP);
  return 0;
}

void printUsage() {
  fprintf(stderr,"Usage: fastasplit [-s/--size] [-f/--file_prefix <prefix>] <fastaFile> <nchunk> <outdir>\n");
}

int strippath(char *fname,char *result) {
  char *chp;
  
  chp=strrchr(fname,'/');
  
  if (chp!=NULL) {
    chp++;
    strcpy(result,chp);
    return (strlen(fname)-strlen(chp));
  } else {
    strcpy(result,fname);
    return 0;
  }
}


int stripext(char *fname,char *result) {
  char *chp;
  int   LenExt;
  
  chp = strrchr(fname,'.');
  
  /* Files which start with a . are treated specially */
  if (chp!=NULL && chp!=fname) {
    LenExt = strlen(fname)-strlen(chp);
    strncpy(result,fname,LenExt);
    return (LenExt);
  } else {
    strcpy(result,fname);
    return 0;
  }
}

