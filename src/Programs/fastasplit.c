#include <stdio.h>
//#include <malloc.h>
#include <stdlib.h>
#include <string.h>

#define MAXSTRLEN 1024*32
#define HEADER_ALLOC_SIZE 1000

int main(int argc, char *argv[]) {
  FILE *InP;
  FILE *OutP = NULL;
  char outdir[MAXSTRLEN];
  char line[MAXSTRLEN];
  long *headerPositions;
  int nHeader = 0;
  long currentPos = 0;
  int i;
  long randPos1;
  long randPos2;
  long tmpPos;
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
  long entryCount;
  int argnum;
  long remainder; 

  if (argc<4) {
    fprintf(stderr,"Usage: fastasplit [-size] <fastaFile> <nchunk> <outdir>\n");
    exit(1);
  } 

  argnum = 1;
  if (!strcmp(argv[1],"-size")) {
    bySize = 1;
    argnum ++;
  }

  
  fileArgNum = argnum++;
  if ((InP = fopen(argv[fileArgNum],"r")) == NULL) {
    fprintf(stderr,"Error: Couldn't open file %s\n",argv[fileArgNum]);
    exit(1);
  }
  strippath(argv[fileArgNum],baseName); 
  stripext(baseName,databaseName); 

  chunkArgNum = argnum++;
  if (!(nChunk = atol(argv[chunkArgNum]))) {
    fprintf(stderr,"Error: Couldn't parse nChunk into a number %s\n",argv[chunkArgNum]);
    exit(1);
  } 


  dirArgNum = argnum++;
  strcpy(outdir,argv[dirArgNum]);


  if ((headerPositions = (long *)calloc(HEADER_ALLOC_SIZE,sizeof(long)))==NULL) {
    fprintf(stderr,"Error: Allocating headerPositions array\n");
    exit(1);
  }

  while (fgets(line,MAXSTRLEN,InP)) {
    if (line[0] == '>') {
      headerPositions[nHeader++] = currentPos;
      if (nHeader && (nHeader%HEADER_ALLOC_SIZE) == 0) {
        if ((headerPositions = (long *)realloc(headerPositions,(nHeader+HEADER_ALLOC_SIZE)*sizeof(long)))==NULL) {
          fprintf(stderr,"Error: Allocating headerPositions array\n");
          exit(1);
        }
      }
    }
    currentPos = ftell(InP);
  }
  if (bySize) { 
    fileSize = currentPos/nChunk;
  
    if (fileSize < 100) {
      fprintf(stderr,
           "Do you really want to split you're file into pieces of %d bytes\n",
           fileSize);
      exit(1);
    }
  } else {
    nPerFile = nHeader/nChunk;
    printf("nPerFile = %d nHeader = %d nChunk = %d\n",nPerFile,nHeader,nChunk);
    if (nPerFile < 1) {
      fprintf(stderr,
           "Do you really want to split you're file into pieces of with less than 1 sequence per file\n");
      exit(1);
    } else {
      remainder = nHeader - nChunk*nPerFile;
      printf ("Remainder = %d\n",nHeader - nChunk*nPerFile);
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

      sprintf(chunkFName,"%s/%s_chunk_%7.7d",
              outdir,databaseName,currentChunk++);

      if (remainder && currentChunk == remainder+1) {
        nPerFile--;
      }

      if ((OutP = fopen(chunkFName,"w")) == NULL) {
        fprintf(stderr,"Error: Couldn't open file %s\n",chunkFName);
        exit(1);
      }
    }
    fseek(InP,headerPositions[i],SEEK_SET);
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

