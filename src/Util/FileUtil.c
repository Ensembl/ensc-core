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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <sys/param.h>
#include <unistd.h>
#include <errno.h>
#include <dirent.h>

#include "FileUtil.h"
#include "StrUtil.h"
#include "Error.h"
#include "Stream.h"
#include "StrUtil.h"

/******************************************************************************/
/* Routine    :                                                               */
/*             FileUtil_countToLine()                                         */
/* Role       :                                                               */
/*             To count the number of lines in a file upto a line beginning   */
/*             with a defined string.                                         */
/* Arguments  :                                                               */
/*             fp      - the file to look in                                  */
/*             NLine   - the number of lines (RETURNED)                       */
/*             string  - the string (if "" then count to end of file)         */
/* Returns    :                                                               */
/*             1 - succes                                                     */
/*             0 - failure                                                    */
/* History    :                                                               */
/*             30/09/93 SMJS  Initial Implementation                          */
/******************************************************************************/
int FileUtil_countToLine(FILE *fp,int *NLine,char *string) {
  char line[MAXSTRLEN];
  int  count = 0;
  int  found = 0;
  int  lenString = strlen(string);

  fgets(line, MAXSTRLEN, fp);
  while (!feof(fp) && !found) {
    count++;
    StrUtil_rmEndSpaces(line);
    if (lenString && !strncmp(line,string,lenString)) {
      found=1;
    }
    fgets(line,MAXSTRLEN,fp);
  }
  *NLine=count;
  return 1;
}

/******************************************************************************/
/* Routine   :                                                                */
/*            FileUtil_exists()                                               */
/* Role      :                                                                */
/*            To check whether a file with a specified name exists            */
/* Arguments :                                                                */
/*            FName    - name of file to check for                            */
/* Returns   :                                                                */
/*            1 - file exists                                                 */
/*            0 - file does not exist (or an error occurred)                  */
/* History   :                                                                */
/*            03/03/98 SMJS  Initial Implementation                           */
/******************************************************************************/
int FileUtil_exists(char *fName) {
  FILE *fp;

  if (fName == NULL) {
    Error_write(ENULLPOINTER,"FileUtil_exists",ERR_SEVERE,"fName is NULL");
    return 0;
  }
  if (strlen(fName) == 0) {
    Error_write(EZEROSTRING,"FileUtil_exists",ERR_SEVERE,"fName");
    return 0;
  }

  if((fp=fopen(fName,"r")) == NULL) {
    return 0;
  }
  fclose(fp);

  return 1;
}

/******************************************************************************/
/* Routine    :                                                               */
/*             FileUtil_getRealPath()                                         */
/* Role       :                                                               */
/*             To get the real path to a file                                 */
/* Arguments  :                                                               */
/*             FNIn  - filename to get path for                               */
/*             Path  - the path to file (RETURNED)                            */
/* Returns    :                                                               */
/*             1 - success                                                    */
/*             0 - failure                                                    */
/* History    :                                                               */
/*             14/02/00 SMJS  Initial Implementation                          */
/* Notes      :                                                               */
/*             realpath() has the nasty habit of changing directories if it   */
/*             encounters an error. That is the reason for all the errno      */
/*             checking and CurDir setting.                                   */
/******************************************************************************/
int FileUtil_getRealPath(char *FNIn, char *Path) {
  char CurDir[MAXPATHLEN];

  if (getcwd(CurDir,MAXPATHLEN)==NULL) {
    Error_write(EGETCWD,"FileUtil_getRealPath",ERR_SEVERE,"errno = %d\n",errno);
    return 0;
  }

  if (errno) {
    Error_write(EGETCWD,"FileUtil_getRealPath",ERR_SEVERE,"errno = %d\n",errno);
    return 0;
  }
   
  if (realpath(FNIn,Path) == NULL) {
    Error_write(EREALPATH,"FileUtil_getRealPath",ERR_SEVERE,
                "realpath returned NULL. File = %s Path = %s.\n",FNIn,Path);
  }

  if (errno && errno!=EINVAL) {
    Error_write(EREALPATH,"FileUtil_getRealPath",ERR_SEVERE,
             "File = %s Path = %s. errno = %d\n",FNIn,Path,errno);
    return 0;
  }

  if (chdir(CurDir)) {
    Error_write(ECHDIR,"FileUtil_getRealPath",ERR_SEVERE,"CurDir = %s. errno = %d\n",CurDir,
             errno);
    return 0;
  }
  return 1;
} 

/******************************************************************************/
/* Routine    :                                                               */
/*             FileUtil_getStrippedLine()                                     */
/* Role       :                                                               */
/*             To get a line from a file, checking that a complete line has   */
/*             been read (by looking for a new line at the end) and then      */
/*             stripping the newline                                          */
/* Arguments  :                                                               */
/*             Line   - the string to fill                                    */
/*             MaxLen - the maximum number of characters to read              */
/*             Fp     - the file to read from.                                */
/* Returns    :                                                               */
/*             1 - success                                                    */
/*             0 - failure                                                    */
/* History    :                                                               */
/*             03/03/98 SMJS  Initial Implementation                          */
/******************************************************************************/
int FileUtil_getStrippedLine(char *line,int maxLen, FILE *Fp) {
  int lenLine;

  line[0]='\0';
  fgets(line,maxLen,Fp);

  lenLine = strlen(line);

  if (lenLine == maxLen-1 && line[lenLine-1] != '\n') {
    Error_write(ESTRTOLONG,"FileUtil_getStrippedLine",ERR_SEVERE,
             "Didn't read a newline\n%s",line);
    return 0;
  }

  StrUtil_stripnewline(line);

  return 1;
}

/******************************************************************************/
/* Routine    :                                                               */
/*             FileUtil_readToChar()                                          */
/* Role       :                                                               */
/*             To read characters upto the first occurence of a particular    */
/*             character or to the end of file if character not found.        */
/* Arguments  :                                                               */
/*             Fp      - the file to get from                                 */
/*             String  - the string (RETURNED)                                */
/*             StrLen  - the maximum length of String                         */
/*             EndChar - the character to use as the termination character    */
/* Returns    :                                                               */
/*             1 - success                                                    */
/*             0 - failure                                                    */
/* History    :                                                               */
/*             11/06/97 SMJS  Initial Implementation                          */
/******************************************************************************/
int FileUtil_readToChar(FILE *Fp,char *string,int strLen,char endChar) {
  int count=0;
  char ch;

  ch=fgetc(Fp);
  while (!feof(Fp) && count<strLen-1 && ch!=endChar) {
    string[count++]=ch;
    ch=getc(Fp);
  }

  if (count == strLen-1) {
    Error_write(ESTRTOLONG,"FileUtil_readToChar",ERR_SEVERE,"strLen = %d",strLen);
    return 0;
  }
  string[count] = '\0';
  return 1;
}

/******************************************************************************/
/* Routine   :                                                                */
/*            FileUtil_close()                                                */
/* Role      :                                                                */
/*            To close a file (and allow debuging)                            */
/* Arguments :                                                                */
/*            fp    - FILE pointer of file to close                           */
/*            fname - name of file to close                                   */
/* Returns   :                                                                */
/*            None                                                            */
/* History   :                                                                */
/*            14/06/93 SMJS  Initial Implementation                           */
/*            10/11/99 SMJS  Added special handling of stdout and stderr      */
/******************************************************************************/
void FileUtil_close(FILE *fp,char *fname) {
#ifdef DBG
  FPrintF(DBGStream,"DBG>> Closing file %s\n",fname);
#endif

  if (strcmp(fname,"stdout") && strcmp(fname,"stderr")) {
    fclose(fp);
  }
}

/******************************************************************************/
/* Routine    :                                                               */
/*             FileUtil_countLines()                                          */
/* Role       :                                                               */
/*             To count the number of lines in a file beginning with a defined*/
/*             string.                                                        */
/* Arguments  :                                                               */
/*             fp      - the file to look in                                  */
/*             string  - the string (if "" then count all lines)              */
/* Returns    :                                                               */
/*             The number of lines beginning with the string                  */
/* History    :                                                               */
/*             14/06/93 SMJS  Initial Implementation                          */
/******************************************************************************/
int FileUtil_countLines(FILE *fp,char *string) {
  int count=0;
  char line[MAXSTRLEN];
  int lenString = strlen(string);

  fgets(line,MAXSTRLEN,fp);
  while (!feof(fp)) {
    if (!strncmp(line,string,lenString)) {
      count ++;
    }
    fgets(line,MAXSTRLEN,fp);
  }
  return count;
}

/******************************************************************************/
/* Routine    :                                                               */
/*             FileUtil_nextPrefixedLine()                                    */
/* Role       :                                                               */
/*             To find the next occurance of a line beginning with a string in*/
/*             a file.                                                        */
/* Arguments  :                                                               */
/*             fp    - the file to look in                                    */
/*             string- the string to look for                                 */
/*             line  - the returned line string                               */
/* Returns    :                                                               */
/*             1 - an occurance of the string was found.                      */
/*             0 - Error eg. EOF                                              */
/* History    :                                                               */
/*             14/06/93 SMJS  Initial Implementation                          */
/******************************************************************************/
int FileUtil_nextPrefixedLine(FILE *fp,char *string,char *line) {
  int lenString = strlen(string);

  fgets(line,MAXSTRLEN,fp);

  while (!feof(fp)) {
    if (!strncmp(line,string,lenString)) {
      return 1;
    }
    fgets(line,MAXSTRLEN,fp);
  }
  return 0;
}
  
/******************************************************************************/
/* Routine   :                                                                */
/*            FileUtil_open()                                                 */
/* Role      :                                                                */
/*            To open a file with error checking                              */
/* Arguments :                                                                */
/*            fname    - name of file to open                                 */
/*            stat     - the mode to open the file in                         */
/*            routine  - where FileUtil_open was called from                  */
/* Returns   :                                                                */
/*            A pointer to the file                                           */
/*            OR null on error                                                */
/* History   :                                                                */
/*            ??/??/?? SMJS  Initial Implementation                           */
/*            10/11/99 SMJS  Added stdout and stderr special handling         */
/******************************************************************************/
FILE *FileUtil_open(char *fname,char *stat,char *routine) {
  FILE *fp;
   
  if (fname == NULL || stat == NULL || routine == NULL) {
    Error_write(ENULLPOINTER,"FileUtil_open",ERR_SEVERE,
             "one or more of fname,stat or routine");
    return NULL;
  }

  if (fname[0] == '\0' || stat[0] == '\0') {
    Error_write(EZEROSTRING,"FileUtil_open",ERR_SEVERE,"one or more of fname or stat");
    return NULL;
  }

  if (!strcmp(fname,"stdout")) {
    fp = stdout;
  } else if (!strcmp(fname,"stderr")) {
    fp = stderr;
  } else if((fp=fopen(fname,stat))==NULL) {
    Error_write(EFOPEN,routine,ERR_SEVERE,"File %s. Mode %s",fname,stat);
    return NULL;
  }
  return fp;
}

/******************************************************************************/
/* Routine    :                                                               */
/*             FileUtil_stripExt()                                            */
/* Role       :                                                               */
/*             To remove the extension from a filename string                 */
/* Arguments  :                                                               */
/*             fname - The filename to remove the extension of                */
/*             result- The stripped filename                                  */
/* Returns    :                                                               */
/*             The number of characters stripped                              */
/* History    :                                                               */
/*             25/03/99 SMJS  Initial Implementation                          */
/******************************************************************************/
int FileUtil_stripExt(char *fname,char *result) {
  char *chp;
  int   LenExt;

  chp = strrchr(fname,'.');

/* Files which start with a . are treated specially */
  if (chp!=NULL && chp!=fname) {
    LenExt = strlen(fname)-strlen(chp);
    StrUtil_strncpy(result,fname,LenExt);
    return (LenExt);
  } else {
    strcpy(result,fname);
    return 0;
  }
}

/******************************************************************************/
/* Routine    :                                                               */
/*             FileUtil_stripPath()                                           */
/* Role       :                                                               */
/*             To remove the path from a filename                             */
/* Arguments  :                                                               */
/*             fname - The filename to remove the path of                     */
/*             result- The stripped filename                                  */
/* Returns    :                                                               */
/*             The number of characters stripped                              */
/* History    :                                                               */
/*             14/06/93 SMJS  Initial Implementation                          */
/******************************************************************************/
int FileUtil_stripPath(char *fname,char *result) {
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
