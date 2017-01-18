/*
 * Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
 * Copyright [2016-2017] EMBL-European Bioinformatics Institute
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

/* 
 * Program to allow moving of a perl or other software installation 
 * 
 * Modifies all occurences of the old location to the new location
 * in the moved files (both text and binary files).
 *
 * For binary file moves the destination path MUST be shorter than
 * the source path.
 *
 * It was written for the otter/apollo/Ensembl website distribution
 * 
 * Written by Steve Searle <searle@sanger.ac.uk>
 */
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <sys/param.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <ctype.h>
#include <errno.h>
#include <string.h>

typedef enum FileTypeEnum {
  NONE,
  EMPTY,
  SCRIPT,
  BINARY,
  TEXT
} FileType;

char *typeStrings[] = {
  "NONE",
  "EMPTY",
  "SCRIPT",
  "BINARY",
  "TEXT"
};

int isDirectory(char *name);
int fileExists(char *name);
void descendAndDo(DIR *dir, char *fromDirName, char *toDirName, char *baseSrcPath,
                  char *baseFromName, char *baseToName, int allowAppend, 
                  char *archiveExtension, int doFunc());
char *copyAndReplace(char *fromBuf,char *toBuf,char *baseSrcName, char *baseToName, char *toName, 
                     int fromLen, int *toLen, FileType type);
FileType fileType(char *fileBuf,int len);
int doFuncInp(char *fromName, char *toName, char *baseSrcName, 
           char *baseToName, int allowAppend, char *archiveExtension);
void printUsageAndDie(void);

int main(int argc, char *argv[]) {
  DIR *fromdir;
  char fullFromName[MAXPATHLEN];
  char fullDestName[MAXPATHLEN];
  char baseSrcPath[MAXPATHLEN];
  char baseFromPath[MAXPATHLEN];
  char baseToPath[MAXPATHLEN];
  struct stat st;
  int  allowAppend = 0;
  int  allowLongerDest = 0;
  char archiveExtension[MAXPATHLEN];

  baseSrcPath[0] = '\0';
  archiveExtension[0] = '\0';

  argv++; argc--;

  while (argc && argv[0][0] == '-') {
    if (strlen(argv[0]) != 2) printUsageAndDie();
    switch (argv[0][1]) {
      case 's':
        argv++; argc--;
        if (!argc || argv[0][0] == '-') printUsageAndDie();
        strcpy(baseSrcPath,argv[0]);
        printf("Source directory set to %s\n",baseSrcPath);
        argv++; argc--;
        break;
      case 'a':
        printf("Allowing append\n");
        allowAppend = 1;
        argv++; argc--;
        break;
      case 'l':
        printf("Allowing destination to be longer than source\n");
        allowLongerDest = 1;
        argv++; argc--;
        break;
      case 'r':
        argv++; argc--;
        if (!argc || argv[0][0] == '-') printUsageAndDie();
        sprintf(archiveExtension,".%s",argv[0]);
        printf("Archive extension set to %s\n",archiveExtension);
        argv++; argc--;
        break;
      default:
        fprintf(stderr,"Error: Unknown option %s\n",argv[0]);
        printUsageAndDie();
        break;
    }
  }

  if (argc != 2) {
    fprintf(stderr,"Error: Didn't end up with two dirs\n");
    printUsageAndDie();
  }
 

  if ((fromdir = opendir(argv[0])) == NULL) {
    fprintf(stderr,"Error: Couldn't open source perl directory %s\n",argv[0]);
    exit(1);
  }

  if ((opendir(argv[1])) != NULL) {
    if (!allowAppend) {
      fprintf(stderr,"Error: Destination perl directory already exists\n");
      exit(1);
    }
  } else if (allowAppend) {
    fprintf(stderr,"Error: -a (append) specified but destination dir doesn't exist\n");
  }

/*
  if (!realpath(argv[0],fullFromName)) {
    fprintf(stderr,"Error: Couldn't get realpath to fromdir\n");
    exit(1);
  }
  if (!realpath(argv[1],fullDestName)) {
    fprintf(stderr,"Error: Couldn't get realpath to destdir\n");
    exit(1);
  }
*/
  strcpy(fullFromName,argv[0]);
  strcpy(fullDestName,argv[1]);


  if (!allowAppend) {
    stat(fullFromName,&st);
    if (mkdir(fullDestName,st.st_mode)) {
      fprintf(stderr, "Error: Failed making directory %s\n",fullDestName);
      exit(1);
    }
  }

  strcpy(baseFromPath,fullFromName);
  strcpy(baseToPath,fullDestName);
  if (baseSrcPath[0]=='\0') {
    strcpy(baseSrcPath,fullFromName);
  }
  if (strlen(baseSrcPath) < strlen(baseToPath) && !allowLongerDest) {
    fprintf(stderr,"Error: Length of to path can't be longer than from path\n");
    exit(1);
  }

  descendAndDo(fromdir,fullFromName,fullDestName,
               baseSrcPath,baseFromPath,baseToPath,
               allowAppend,archiveExtension,doFuncInp);
}

void printUsageAndDie() {
  fprintf(stderr, "Usage: moveperl [-l] [-a] [-s <srcdir>] [-r <libsuffix>] <fromdir> <todir>\n");
  exit(1);
}

void descendAndDo(DIR *dir, char *fromDirName, char *toDirName, char *baseSrcName,
                  char *baseFromName, char *baseToName, int allowAppend, 
                  char *archiveExtension, int doFunc()) {
  struct dirent *dp;
  char fromName[MAXPATHLEN];
  char toName[MAXPATHLEN];
  struct stat st;

  while ((dp = readdir(dir)) != NULL) {
    sprintf(fromName,"%s/%s",fromDirName,dp->d_name);
    sprintf(toName,"%s/%s",toDirName,dp->d_name);
    
    if (strcmp(dp->d_name,".") && strcmp(dp->d_name,"..")) {
      stat(fromName,&st);
      if (st.st_mode & S_IFDIR) {
        DIR *childdir;
        //printf("  is a directory\n");
        if ((childdir = opendir(fromName)) == NULL) {
          fprintf(stderr,"Error: Couldn't open directory %s\n",fromName);
          exit(1);
        }
        if (!doFunc(fromName,toName,baseSrcName,baseToName,
                    allowAppend,archiveExtension)) {
          fprintf(stderr,"Error: Failed executing doFunc\n");
          exit(1);
        }
/* Recurse here */
        descendAndDo(childdir,fromName,toName,baseSrcName,baseFromName,
                     baseToName,allowAppend,archiveExtension,doFunc);
      } else {
        if (!doFunc(fromName,toName,baseSrcName,baseToName,
                    allowAppend,archiveExtension)) {
          fprintf(stderr,"Error: Failed executing doFunc\n");
          exit(1);
        }
      }
    } 
  }
  (void)closedir(dir);
}

int doFuncInp(char *fromName, char *toName, char *baseSrcName, 
           char *baseToName, int allowAppend, char *archiveExtension) {
  FILE *fp;
  char *fileBuf;
  char *toBuf;
  int ch;
  char *chP;
  struct stat st;
  int count=0;
  FileType type;
  int toLen;
  int toNameLen = strlen(toName);
  int archExtLen = strlen(archiveExtension);
  
  //printf("%s %s %s %s\n",fromName,toName,baseSrcName,baseToName);

  stat(fromName,&st);
  if (st.st_mode & S_IFDIR) {
    if (!fileExists(toName)) {
      if (mkdir(toName,st.st_mode)) {
        fprintf(stderr,"Error: Failed making directory %s\n",toName);
        exit(1);
      }
    } else if (!isDirectory(toName)) {
      fprintf(stderr,"Error: %s is not a directory, but should be\n",toName);
      exit(1);
    }
  } else {
    if (fileExists(toName)) {
      fprintf(stderr,"Error: File %s already exists and I'm not going to overwrite it!",toName);
      exit(1);
    }
    if ((fileBuf = calloc(st.st_size+1,sizeof(char))) == NULL) {
      fprintf(stderr,"Error: Failed allocating space for file buffer of size %ld\n", st.st_size+1);
      exit(1);
    }
/* Hack */
    if ((toBuf = calloc(st.st_size+1000000,sizeof(char))) == NULL) {
      fprintf(stderr,"Error: Failed allocating space for file buffer\n");
      exit(1);
    }
  
    if ((fp = fopen(fromName,"r")) == NULL) {
      fprintf(stderr,"Error: Failed opening from perl file %s\n",fromName);
      exit(1);
    }
    fread(fileBuf,st.st_size,1,fp); 
/*
    chP = fileBuf;
    
    while ((ch = getc(fp)) != EOF) {
      count++;
      *chP = ch;
      chP++;
    }
*/
    fclose(fp);

    if ((fp = fopen(toName,"w")) == NULL) {
      fprintf(stderr,"Error: Failed opening to perl file %s\n",toName);
      exit(1);
    }
    type = fileType(fileBuf,st.st_size);
    //printf("Type = %d (%s)\n",type,typeStrings[type]);
    copyAndReplace(fileBuf,toBuf,baseSrcName,baseToName,toName,st.st_size,&toLen,type);
    fwrite(toBuf,toLen,1,fp); 
    fclose(fp);
    chmod(toName,st.st_mode);

    
    if (archExtLen && toNameLen > archExtLen) {

      //printf("toNameLen = %d archExtLen = %d toName bit = %s\n",
      //       toNameLen,archExtLen,&(toName[toNameLen-archExtLen])); 

      if (!strcmp(&(toName[toNameLen-archExtLen]),archiveExtension)) {
        char comStr[MAXPATHLEN + 1024];
        sprintf(comStr,"ranlib %s",toName);
        printf("%s\n",comStr);
        system(comStr);
      }
    }

    free(fileBuf);
    free(toBuf);
  }
  return 1;
}


int fileExists(char *name) {
  struct stat st;

  if (stat(name,&st)) {
    if (errno == ENOENT) {
      return 0;
    } else {
      fprintf(stderr,"Warning: Stat of %s failed for reason other than not found\n",name);
      perror(NULL);
      return 0;
    }
  }
  return 1;
}

int isDirectory(char *name) {
  struct stat st;

  if (fileExists(name)) {
    stat(name,&st);
    if (st.st_mode & S_IFDIR) {
      return 1;
    }
  }
  return 0;
}

char *copyAndReplace(char *fromBuf,char *toBuf,char *baseSrcName, char *baseToName, char *toName, 
                     int fromLen, int *toLen, FileType type) {
  char *fromChP;
  char *toChP;
  int nToCopy = 0;
  int lenBaseSrcName = strlen(baseSrcName);
  int lenBaseToName = strlen(baseToName);
  int i;
  int printDone = 0;
  int lenPad = 0;

  
  
  fromChP = fromBuf;
  toChP = toBuf;


  //printf("base source name = %s\n",baseSrcName);
  while ((fromChP-fromBuf) < fromLen) {
    if (*fromChP == baseSrcName[0] && 
        !strncmp(baseSrcName,fromChP,lenBaseSrcName)) {
      strcpy(toChP,baseToName);

    
      if (!printDone) {
        printf("Modifying %s type %s\n",toName,typeStrings[type]);
        printDone = 1;
      }

      fromChP  += lenBaseSrcName;
      toChP    += lenBaseToName;
      lenPad   = lenBaseSrcName-lenBaseToName;
      while (!iscntrl(*fromChP) && !isspace(*fromChP) && !(*fromChP==':' && type == BINARY)) {
        //printf("%c",*fromChP);
        *toChP = *fromChP;
        toChP++;
        fromChP++;
      }
      if (type == BINARY) {
        // Handling ':' separated lists in binary files - need to loop through
        // all list elements, keeping track of how much padding we need to add
        // at the end of the list
        while (*fromChP==':') {
          *toChP = *fromChP;
          toChP++;
          fromChP++;

          if (*fromChP == baseSrcName[0] && 
             !strncmp(baseSrcName,fromChP,lenBaseSrcName)) {
            lenPad   += lenBaseSrcName-lenBaseToName;
            strcpy(toChP,baseToName);
            fromChP  += lenBaseSrcName;
            toChP    += lenBaseToName;
          }
          while (!iscntrl(*fromChP) && !isspace(*fromChP) && !(*fromChP==':' && type == BINARY)) {
            //printf("%c",*fromChP);
            *toChP = *fromChP;
            toChP++;
            fromChP++;
          }
        }
        //printf("padding with %d chars\n",lenPad);
        for (i=0;i<lenPad;i++) {
          (*toChP)='\0';
          toChP++;
        }
        lenPad = 0;
      }
      //printf("\n");
    } else {
      *toChP = *fromChP;
      fromChP++;
      toChP++;
    }
  } 
  
  *toLen = toChP-toBuf;
  if (type == BINARY && fromLen!=*toLen) {
    fprintf(stderr, "Error: BINARY file ended up different length after modification\n");
    fprintf(stderr, "From len = %d,  To len = %d\n",fromLen, *toLen);
    exit(1);
  }
  //printf("From len = %d,  To len = %d\n",fromLen, *toLen);
  return toBuf;
}

FileType fileType(char *fileBuf,int len) {
  char *chP;
  int i;
  int nNewline = 0;
  int nCntrl = 0;
  int nToGet=2000;

  if (!len) return EMPTY;

  if (!strncmp(fileBuf,"#!",2)) {
    return SCRIPT;  
  }

  if (len < nToGet) {
    nToGet = len;
  }
  for (i=0;i<nToGet;i++) {
    if (fileBuf[i] == '\n') {
      nNewline++;
    } else if (iscntrl(fileBuf[i]) && fileBuf[i]!='\t') {
      nCntrl++;
    }
  }
  //printf("nNewline = %d nCntrl = %d nToGet = %d\n",nNewline,nCntrl,nToGet);
  if (!nNewline) return BINARY;

  if ((float)nCntrl/(float)nToGet > 0.15) return BINARY;

  if ((float)nCntrl/(float)nToGet < 0.01) return TEXT;
  if ((float)nToGet/(float)nNewline < 81) return TEXT;

  return BINARY;
}
