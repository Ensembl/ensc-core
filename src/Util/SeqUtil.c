#include "SeqUtil.h"
#include "StrUtil.h"
#include "Error.h"
#include "Stream.h"
#include "FileUtil.h"

char *SeqUtil_reverseComplement(char *seqStr, int lenSeqStr) {
  int i;

  StrUtil_reverseString(seqStr,lenSeqStr);

  for (i=0;i<lenSeqStr;i++) {
    switch (seqStr[i]) {
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
        fprintf(stderr,"ERROR: Failed reverse complementing\n");
        return NULL;
    }
  }
  
  return seqStr;
}

char *SeqUtil_addNs(char *seq, int length) {
// NIY Make fast
  int seqlen = strlen(seq);
  char *chP;
  int i;

  if ((seq = (char *)realloc(seq,(strlen(seq)+1+length)*sizeof(char))) == NULL) {
    fprintf(stderr, "ERROR: Failed reallocating sequence string\n");
    exit(1);
  }
  chP = &(seq[seqlen]);
  for (i=seqlen; i<seqlen+length; i++) {
    *chP = 'N';
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
  char  TmpName[MAXSTRLEN];
  char *headerStr = "Symbol   3-letter    Codons";

  if ((FpIn=FileUtil_open(fName,"r","ReadTransTab"))==NULL) {
    Error_trace("ReadTransTab",NULL);
    return 0;
  }

  fgets(Line,MAXSTRLEN,FpIn);

  while (strncmp(Line,headerStr,sizeof(headerStr)) && !feof(FpIn)) {
    fgets(Line,MAXSTRLEN,FpIn);
  }

  fgets(Line,MAXSTRLEN,FpIn);

  while(!feof(FpIn)) {
    fgets(Line,MAXSTRLEN,FpIn);
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
        if (!NucToInt(Token[i],&(Inds[i]))) {
          Error_trace("ReadTransTab",NULL);
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
int NucToInt(char Nuc, int *NucInt) {
   switch (Nuc)
   {
      case 'A':
      case 'a':
         *NucInt = 0;
         break;
      case 'C':
      case 'c':
         *NucInt = 1;
         break;
      case 'G':
      case 'g':
         *NucInt = 2;
         break;
      case 'T':
      case 't':
         *NucInt = 3;
         break;
      default:
         Error_write(EFALLTHRU,"NucToInt",ERR_SEVERE,
                     "Nucleotide = %c (not A,C,G or T)",Nuc);
         return 0;
   }
   return 1;
}

