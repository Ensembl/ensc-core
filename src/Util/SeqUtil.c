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

#define __SEQUTIL_MAIN__
#include "SeqUtil.h"
#undef __SEQUTIL_MAIN__
#include "StrUtil.h"
#include "Error.h"
#include "Stream.h"
#include "FileUtil.h"

#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>

char *SeqUtil_reverseComplement(char *seqStr, int lenSeqStr) {
  int i;

  StrUtil_reverseString(seqStr,lenSeqStr);

  for (i=0;i<lenSeqStr;i++) {
    switch (toupper(seqStr[i])) {
      case 'A':
        seqStr[i]='T';
        break;
      case 'C':
        seqStr[i]='G';
        break;
      case 'G':
        seqStr[i]='C';
        break;
      case 'T':
        seqStr[i]='A';
        break;
      case 'N':
        seqStr[i]='N';
        break;
      case 'R':
        seqStr[i]='Y';
        break;
      case 'Y':
        seqStr[i]='R';
        break;
      case 'M':
        seqStr[i]='K';
        break;
      case 'K':
        seqStr[i]='M';
        break;
      case 'S':
        seqStr[i]='W';
        break;
      case 'W':
        seqStr[i]='S';
        break;
      case 'H':
        seqStr[i]='D';
        break;
      case 'D':
        seqStr[i]='H';
        break;
      case 'B':
        seqStr[i]='V';
        break;
      case 'V':
        seqStr[i]='B';
        break;
      default:
        fprintf(stderr,"ERROR: Failed reverse complementing char = %c\n", seqStr[i]);
        return NULL;
    }
  }
  
  return seqStr;
}

char *SeqUtil_addGaps(char *seq, int length) {
  return SeqUtil_addRes(seq,length,'-');
}

char *SeqUtil_addNs(char *seq, int length) {
  return SeqUtil_addRes(seq,length,'N');
}

char *SeqUtil_addRes(char *seq, int length, char res) {
// NIY Make fast
  int seqlen = strlen(seq);
  char *chP;
  int i;

  //printf("Seqlen = %d length = %d\n",seqlen,length);
  if (length == -1) return seq;

  if ((seq = (char *)realloc(seq,(seqlen+1+length)*sizeof(char))) == NULL) {
    fprintf(stderr, "ERROR: Failed reallocating sequence string for length %d pointer %p\n",seqlen+1+length, seq);
    exit(1);
  }
  chP = &(seq[seqlen]);
  for (i=seqlen; i<seqlen+length; i++) {
    *chP = res;
    chP++;
  }
  *chP = '\0';
  
  return seq;
}


/******************************************************************************/
/* Routine    :                                                               */
/*             ReadTransTab()                                                 */
/* Role       :                                                               */
/* Arguments  :                                                               */
/* Returns    :                                                               */
/* History    :                                                               */
/******************************************************************************/
int SeqUtil_readTransTab(char *fName, char TransTab[4][4][4]) {
  FILE *FpIn;
  char  Line[MAXSTRLEN];
  char  Token[MAXSTRLEN];
  int   i;
  int   j;
  int   k;
  int   Inds[3];
  char  aa;
  char *ChP;
  char *headerStr = "Symbol   3-letter    Codons";

  if ((FpIn=FileUtil_open(fName,"r","ReadTransTab"))==NULL) {
    Error_trace("ReadTransTab",NULL);
    return 0;
  }

  if  (!fgets(Line,MAXSTRLEN,FpIn)) {
    Error_write(EFGETSNULL,"ReadTransTab", ERR_SEVERE, "reading translation table file %s\n", fName);
    return 0;
  }
  

  while (strncmp(Line,headerStr,strlen(headerStr)) && !feof(FpIn)) {
    if  (!fgets(Line,MAXSTRLEN,FpIn)) {
      Error_write(EFGETSNULL,"ReadTransTab", ERR_SEVERE, "reading translation table file %s\n", fName);
      return 0;
    }
  }

  if  (!fgets(Line,MAXSTRLEN,FpIn)) {
    Error_write(EFGETSNULL,"ReadTransTab", ERR_SEVERE, "reading translation table file %s\n", fName);
    return 0;
  }

  while(!feof(FpIn)) {
    if  (!fgets(Line,MAXSTRLEN,FpIn)) {
      Error_write(EFGETSNULL,"ReadTransTab", ERR_SEVERE, "reading translation table file %s\n", fName);
      return 0;
    }
    ChP = Line;
    if (!StrUtil_gettok(Token,&ChP,ChP,MAXSTRLEN)) {
      Error_trace("ReadTransTab",NULL);
      return 0;
    }
    if (strlen(Token)!=1) {
      Error_write(ESTRTOLONG,"ReadTransTab",ERR_SEVERE,
               "String should be aa single letter code = %s\n",Token);
      return 0;
    }
    aa = Token[0];
    if (!StrUtil_gettok(Token,&ChP,ChP,MAXSTRLEN)) {
      Error_trace("ReadTransTab",NULL);
      return 0;
    }
    if (strlen(Token)!=3) {
      Error_write(ESTRTOLONG,"ReadTransTab",ERR_SEVERE,
               "String should be aa three letter code = %s\n",Token);
      return 0;
    }
    if (!StrUtil_gettok(Token,&ChP,ChP,MAXSTRLEN)) {
      Error_trace("ReadTransTab",NULL);
      return 0;
    }
    while (strcmp(Token,"!") && strcmp(Token,"...")) {
      if (strlen(Token)!=3) {
        Error_write(ESTRTOLONG,"ReadTransTab",ERR_SEVERE,
                 "String should be a three letter nucleotide = %s\n",Token);
        return 0;
      }
      for (i=0;i<3;i++) {
        if ((Inds[i] = nucToIntArray[Token[i]]) == -1) {
          Error_trace("ReadTransTab",Token);
          SeqUtil_printConvTable(nucToIntArray);
          return 0;
        }
      }
      TransTab[Inds[0]][Inds[1]][Inds[2]] = aa;
      if (!StrUtil_gettok(Token,&ChP,ChP,MAXSTRLEN)) {
        Error_trace("ReadTransTab",NULL);
        return 0;
      }
    }
  }
  for (i=0;i<4;i++) {
    for (j=0;j<4;j++) {
      Stream_fprintf(OutStream,"{");
      for (k=0;k<4;k++) {
        Stream_fprintf(OutStream,"%c ",TransTab[i][j][k]);
      }
      Stream_fprintf(OutStream,"} ");
    }
    Stream_fprintf(OutStream,"\n");
  }
  return 1;
}

void SeqUtil_printConvTable(int *convTable) {
  int i;
  for (i=0; i<256; i++) {
    printf("%d %c %d\n",i, i,convTable[i]);
  }
}

void SeqUtil_writeFasta(FILE *fp, char *header, char *seq, int lineLen) {
  fprintf(fp,">%s\n",header);

  if (!lineLen) lineLen = 60;
  char *chP = seq;
  int linePos =0;
  while (*chP != '\0') {
    fputc(*chP++, fp);
    linePos++;
    if (linePos==lineLen) {
      fputc('\n', fp);
      linePos = 0;
    }
  }
  if (linePos) {
    fputc('\n',fp);
  }
}
