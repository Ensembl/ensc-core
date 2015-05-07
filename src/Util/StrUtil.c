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

#include <ctype.h>
#include <string.h>
#include <stdlib.h>

#define __STRUTIL_MAIN__
#include "StrUtil.h"
#undef __STRUTIL_MAIN__
#include "Error.h"

char *StrUtil_appendString(char *to, char *from) {
  int lenTo = strlen(to);
  int lenFrom = strlen(from);

  //fprintf(stderr,"lenTo = %d lenFrom = %d\n",lenTo, lenFrom);
  if ((to = (char *)realloc(to,lenTo+lenFrom+1))== NULL) {
    fprintf(stderr,"ERROR: Failed reallocating string\n");
    exit(1);
  }

  memcpy(&to[lenTo], from, lenFrom);
//  strcat(to,from);
  to[lenTo+lenFrom] = '\0';

  return to;
}

char *StrUtil_appendNString(char *to, char *from, int len) {
  int lenTo = strlen(to);
  int lenFrom = strlen(from);

  lenFrom = (len > lenFrom) ? lenFrom : len;

  if ((to = (char *)realloc(to,lenTo+lenFrom+1))== NULL) {
    fprintf(stderr,"ERROR: Failed reallocating string\n");
    exit(1);
  }

  strncat(to,from,lenFrom);
  to[lenTo+lenFrom] = '\0';

  return to;
}

void StrUtil_reverseString(char *string, int len) {
  int i;
  char *chPUp,*chPDown;
  char tempCh;
  
  chPUp=&(string[0]);
  chPDown=&(string[len-1]);
  for(i=0;i<len/2;i++) {
    tempCh   = *chPDown;
    *chPDown = *chPUp;
    *chPUp   = tempCh;
    chPUp++;
    chPDown--;
  }
  return;
}

/******************************************************************************/
/* Routine   :                                                                */
/*            StrUtil_strupr()                                                */
/* Role      :                                                                */
/*            To uppercase a string                                           */
/* Arguments :                                                                */
/*            string - string to uppercase                                    */
/* Returns   :                                                                */
/*            Length of string                                                */
/******************************************************************************/
int StrUtil_strupr(char *string) {
  int len=0;
  char *chP = string;

  while (*chP != '\0') {
    *chP = toupper(*chP);
    //len++;
    chP++;
  }
  return len;
}

/******************************************************************************/
/* Routine    :                                                               */
/*             StrUtil_copyString()                                           */
/* Role       :                                                               */
/*             To copy a string from another string starting at a specified   */
/*             position. If the string to be copied is not long enough then   */
/*             allocate a 1 element string and insert string terminator       */
/* Arguments  :                                                               */
/*             to - char pointer to copy to  (RETURNED ALLOCATED)             */
/*             from - string to copy                                          */
/*             Start - the starting position in the from string to copy       */
/* Returns    :                                                               */
/*             1 - success                                                    */
/*             0 - failure                                                    */
/* History    :                                                               */
/*             21/06/94 SMJS  Initial Implementation                          */
/******************************************************************************/
char *StrUtil_copyString(char **to, char *from,int Start) {
  char *cp;
  int   Len;

  if ((Len=strlen(from))-1 < Start) {
    if ((*to = (char *)calloc(1,sizeof(char)))==NULL) {
      Error_write(EMALLERR,"StrUtil_copyString",ERR_SEVERE,"to");
      return NULL;
    }
    (*to)[0]='\0';
  } else {
    Len-=Start;
    cp = &from[Start];
    if ((*to = (char *)calloc(Len+3,sizeof(char)))==NULL) {
      Error_write(EMALLERR,"StrUtil_copyString",ERR_SEVERE,"to");
      return NULL;
    }
    strcpy(*to,cp);
  }
  return *to;
}

/******************************************************************************/
/* Routine   :                                                                */
/*            StrUtil_strlwr()                                                */
/* Role      :                                                                */
/*            To lowercase a string                                           */
/* Arguments :                                                                */
/*            string - string to lowercase                                    */
/* Returns   :                                                                */
/*            Length of string                                                */
/******************************************************************************/
int StrUtil_strlwr(char *string) {
  int len=0;
  char *chP = string;

  while (*chP != '\0') {
    *chP = tolower(*chP);
    len++;
    chP++;
  }
  return len;
}

/******************************************************************************/
/* Routine    :                                                               */
/*             StrUtil_rmTrailSpace()                                         */
/* Role       :                                                               */
/*             Remove trailing spaces from a string (including tabs)          */
/* Arguments  :                                                               */
/*             Str - the string                                               */
/* Returns    :                                                               */
/*             The length of the edited string                                */
/* History    :                                                               */
/*             05/06/93                                                       */
/******************************************************************************/
int StrUtil_rmTrailSpace(char *str) {
  int pos=0;

/* Initialise pos to end of string */
  pos = strlen(str)-1;

  while (pos!=(-1) && (str[pos]==' ' || str[pos]=='\t' || str[pos]=='\n')) {
    pos--;
  }

  str[pos+1]='\0';
  return (pos+1);
}

/******************************************************************************/
/* Routine    :                                                               */
/*             StrUtil_rmLeadSpace()                                          */
/* Role       :                                                               */
/*             To remove any leading spaces from a string (including tabs)    */
/* Arguments  :                                                               */
/*             str - the string                                               */
/* Returns    :                                                               */
/*             The length of the new string                                   */
/* History    :                                                               */
/*             05/06/93 SMJS  Initial Implementation                          */
/*             24/11/95 SMJS  Added check for TrueStart==0                    */
/******************************************************************************/
int StrUtil_rmLeadSpace(char *str) {
  int trueStart=0;
  int len;

  if ((len=strlen(str))) {
    while (trueStart!=len &&
           (str[trueStart]==' ' ||
           str[trueStart]=='\t')) {
      trueStart++;
    }
    if (trueStart!=len && trueStart) {
      memmove(str,&str[trueStart],len-trueStart+1);
    } else if (trueStart) {
      str[0]='\0';
    }
  }

  return (len-trueStart);
}

/******************************************************************************/
/* Routine    :                                                               */
/*             StrUtil_rmEndSpaces()                                          */
/* Role       :                                                               */
/*             To remove leading and trailing spaces from a string            */
/* Arguments  :                                                               */
/*             str - the string to remove spaces from.                        */
/* Returns    :                                                               */
/*             The length of the new string.                                  */
/* History    :                                                               */
/*             04/06/93 SMJS  Initial implementation                          */
/******************************************************************************/
int StrUtil_rmEndSpaces(char *str) {
  StrUtil_rmLeadSpace(str);

  return StrUtil_rmTrailSpace(str);
}

/******************************************************************************/
/* Routine    :                                                               */
/*             StrUtil_strncpy()                                              */
/* Role       :                                                               */
/*             A strncpy which always appends a NULL to the string            */
/* Arguments  :                                                               */
/*             s1    - string to copy to                                      */
/*             s2    - string to copy from                                    */
/*             len   - the number of characters to copy                       */
/* Returns    :                                                               */
/*             a pointer to the new string OR                                 */
/*             NULL on error                                                  */
/* History    :                                                               */
/*             10/06/93 SMJS  Initial Implementation for parse                */
/******************************************************************************/
char *StrUtil_strncpy(char *s1,char *s2,int len) {
  if (len == 0) {
    s1[0]='\0';
    return s1;
  }

  if (s1 == NULL || s2 == NULL) {
    Error_write(ENULLPOINTER,"StrUtil_strncpy",ERR_SEVERE,"s1 or s2");
    return NULL;
  }

  if (strncpy(s1,s2,len)==NULL) {
    Error_write(ESTRNCPYNULL,"StrUtil_strncpy",ERR_SEVERE,"");
    return NULL;
  }

  s1[len]='\0';

  return s1;
}


/******************************************************************************/
/* Routine    :                                                               */
/*             StrUtil_fastCopyString                                         */
/* Role       :                                                               */
/*             To copy a string from another string allocating space for it   */
/*             when the length of the source string is known in advance.      */
/*             This routine does no checks on the strings (use                */
/*             StrUtil_copyNString for that). This routine aims to be         */
/*             efficient.                                                     */
/* Arguments  :                                                               */
/*             to     - char pointer to copy to (RETURNED ALLOCATED)          */
/*             from   - string to copy                                        */
/*             length - the (maximum) number of characters to copy            */
/* Returns    :                                                               */
/*             pointer to string OR                                           */
/*             NULL on failure                                                */
/* History    :                                                               */
/*             11/11/99 SMJS  Initial Implementation                          */
/******************************************************************************/
char *StrUtil_fastCopyString(char **to, char *from,int length) {

  if ((*to = (char *)calloc(length+1,sizeof(char)))==NULL) {
    Error_write(ECALLERR,"StrUtil_copyLString",ERR_SEVERE,"to");
    return NULL;
  }
  memcpy(*to,from,length+1);

  return *to;
}

/******************************************************************************/
/* Routine    :                                                               */
/*             StrUtil_copyNString()                                          */
/* Role       :                                                               */
/*             To copy a string from another string starting at a specified   */
/*             position. If the string to be copied is not long enough then   */
/*             allocate a 1 element string and insert string terminator       */
/* Arguments  :                                                               */
/*             to     - char pointer to copy to (RETURNED ALLOCATED)          */
/*             from   - string to copy                                        */
/*             Start  - the starting position in the from string to copy      */
/*             length - the (maximum) number of characters to copy            */
/* Returns    :                                                               */
/*             pointer to string OR                                           */
/*             NULL on failure                                                */
/* History    :                                                               */
/*             30/08/98 SMJS  Initial Implementation                          */
/******************************************************************************/
char *StrUtil_copyNString(char **to, char *from,int Start,int length) {
  int lenFrom = strlen(from);

  if (lenFrom-1 < Start) {

    if ((*to = (char *)calloc(1,sizeof(char)))==NULL) {
      Error_write(ECALLERR,"StrUtil_copyNString",ERR_SEVERE,"to");
      return NULL;
    }
    (*to)[0]='\0';

  } else {

    if ((*to = (char *)calloc(length+1,sizeof(char)))==NULL) {
      Error_write(ECALLERR,"StrUtil_copyNString",ERR_SEVERE,"*to");
      return NULL;
    }

    StrUtil_strncpy(*to,&(from[Start]),length);

    if ((*to = (char *)realloc(*to, (strlen(*to)+1)*sizeof(char)))==NULL) {
      Error_write(ECALLERR,"StrUtil_copyNString",ERR_SEVERE,"to");
      return NULL;
    }

  }
  return *to;
}

/******************************************************************************/
/* Routine    :                                                               */
/*             StrUtil_copyNTString()                                         */
/* Role       :                                                               */
/*             To copy a string from another string starting at a specified   */
/*             position. If the string to be copied is not long enough then   */
/*             allocate a 1 element string and insert string terminator       */
/* Arguments  :                                                               */
/*             to     - char pointer to copy to (RETURNED ALLOCATED)          */
/*             from   - string to copy                                        */
/*             Start  - the starting position in the from string to copy      */
/*             length - the (maximum) number of characters to copy            */
/* Returns    :                                                               */
/*             pointer to string OR                                           */
/*             NULL on failure                                                */
/* History    :                                                               */
/*             30/08/98 SMJS  Initial Implementation                          */
/******************************************************************************/
char *StrUtil_copyNTString(char **to, char *from,int Start,int length) {
  int lenTo;
  int lenFrom = strlen(from);

  if (lenFrom-1 < Start) {

    if ((*to = (char *)malloc(sizeof(char)))==NULL) {
      Error_write(ECALLERR,"StrUtil_copyNTSring",ERR_SEVERE,"to");
      return NULL;
    }
    (*to)[0]='\0';

  } else {
    if ((*to = (char *)calloc(length+1,sizeof(char)))==NULL) {
      Error_write(ECALLERR,"StrUtil_copyNTString",ERR_SEVERE,"*to");
      return NULL;
    }

    StrUtil_strncpy(*to, &(from[Start]),length);
    lenTo = StrUtil_rmEndSpaces(*to);

    if ((*to = (char *)realloc(*to,(lenTo+1)*sizeof(char)))==NULL) {
      Error_write(EREALLERR,"StrUtil_copyNTtring",ERR_SEVERE,"to");
      return NULL;
    }
  }
  return *to;
}

/******************************************************************************/
/* Routine    :                                                               */
/*             StrUtil_copyTString()                                          */
/* Role       :                                                               */
/*             To copy a string from another string starting at a specified   */
/*             position. If the string to be copied is not long enough then   */
/*             allocate a 1 element string and insert string terminator       */
/* Arguments  :                                                               */
/*             to     - char pointer to copy to (RETURNED ALLOCATED)          */
/*             from   - string to copy                                        */
/*             Start  - the starting position in the from string to copy      */
/* Returns    :                                                               */
/*             pointer to string OR                                           */
/*             NULL on failure                                                */
/* History    :                                                               */
/*             30/08/98 SMJS  Initial Implementation                          */
/******************************************************************************/
char *StrUtil_copyTString(char **to, char *from,int Start) {
  int lenTo;
  int lenFrom = strlen(from);

  if (lenFrom-1 < Start) {
    if ((*to = (char *)calloc(1,sizeof(char)))==NULL) {
      Error_write(ECALLERR,"StrUtil_copyTString",ERR_SEVERE,"to");
      return NULL;
    }
    (*to)[0]='\0';

  } else {

    if ((*to = (char *)calloc(lenFrom+1,sizeof(char)))==NULL) {
      Error_write(ECALLERR,"StrUtil_copyTString",ERR_SEVERE,"*to");
      return NULL;
    }

    strcpy(*to,&(from[Start]));
    lenTo = StrUtil_rmEndSpaces(*to);

    if ((*to = (char *)realloc(*to,(lenTo+1)*sizeof(char)))==NULL) {
      Error_write(EREALLERR,"StrUtil_copyTString",ERR_SEVERE,"to");
      return NULL;
    }
  }
  return *to;
}

/******************************************************************************/
/* Routine    :                                                               */
/*             StrUtil_stripnewline()                                         */
/* Role       :                                                               */
/*             strip newlines from a string                                   */
/* Arguments  :                                                               */
/*             str - the string                                               */
/* Returns    :                                                               */
/*             1 - success                                                    */
/*             0 - failure                                                    */
/* History    :                                                               */
/*             ??/??/?? SMJS  Initial Implementation                          */
/*             10/02/98 SMJS  Changed so that zero length string is not an    */
/*                            error.                                          */
/******************************************************************************/
int StrUtil_stripnewline(char *str) {
  char *ChP;
  int j=0;
  int count=0;


  for (ChP=str;*ChP!='\0';ChP++) {
    if (*ChP!='\n') {
      str[j]= *ChP;
      j++;
    } else {
      count++;
    }
  }
  str[j]='\0';

  return 1;
}

char *StrUtil_strReplChr(char *str, char fromCh, char toCh) {
  char *ChP;

  for (ChP=str;*ChP!='\0';ChP++) {
    if (*ChP==fromCh) {
      *ChP = toCh;
    }
  }
  return str;
}

char *StrUtil_strReplChrs(char *str, char *fromChs, char *toChs) {
  char *chP;
  int nFrom = strlen(fromChs);
  int nTo   = strlen(toChs);
  int i;

  if (nFrom != nTo) {
    fprintf(stderr, "Error: Different numbers of from and to chars in strReplChrs\n");
    exit(1);
  }

  for (i=0; i<nFrom; i++) {
    charMap[fromChs[i]] = toChs[i];
  }

  for (chP=str;*chP!='\0';chP++) {
    *chP = charMap[*chP];
  }

// Revert up the charMap
  for (i=0; i<nFrom; i++) {
    charMap[fromChs[i]] = fromChs[i];
  }

  return str;
}

/******************************************************************************/
/* Routine   :                                                                */
/*            gettok()                                                        */
/* Role      :                                                                */
/*            To get the first white space seperated string from a supplied   */
/*            char array                                                      */
/* Arguments :                                                                */
/*            str    - the extracted token (RETURNED)                         */
/*            out_pp - a pointer to a string pointer to return the position   */
/*                     after the newly read token (RETURNED)                  */
/*            in_p   - the input string                                       */
/*            MaxLen - maximum length for token                               */
/* Returns   :                                                                */
/*            1 - success                                                     */
/*            0 - failure                                                     */
/* Notes     :                                                                */
/*            Tidied up file created 4/93 for human                           */
/* History   :                                                                */
/*            07/06/93 SMJS  Version created for my parser                    */
/*            11/06/93 SMJS  Changed so that it doesn't give an error if no   */
/*                           argument read.                                   */
/*            28/01/94 SMJS  Changed so that can deal with quoted strings     */
/*            31/07/98 SMJS  Added maximum length argument to gettok()        */
/*                           instead of using fixed maximum length            */
/******************************************************************************/
int StrUtil_gettok(char *str,char **out_pp,char *in_p,int MaxLen) {
  int i;
  int NQuote=0;
  int LenInStr;

  if (in_p==NULL) {
    Error_write(ENULLPOINTER,"gettok",ERR_SEVERE,"in_p");
    str[0]='\0';
    return 0;
  }

  i=0;
  LenInStr = strlen(in_p);
  while ((in_p[i]==' ' || in_p[i]=='\t' || in_p[i]=='\n') && i<LenInStr) {
    i++;
  }

  if (i==strlen(in_p)) {
    *out_pp=(&in_p[i]); /* out_pp contains a pointer to the position in */
                        /* the string */
    str[0]='\0';
    return 1;
  }

  if (in_p[i]!='\"') {
    in_p=(&in_p[i]);
    LenInStr = strlen(in_p);
    i=0;
    while (in_p[i]!=' ' && i<MaxLen && in_p[i]!='\0' && in_p[i]!='\n' &&
           in_p[i]!='\t' && i<LenInStr) {
      str[i]=in_p[i];
      i++;
    }
    *out_pp=(&in_p[i]);/* out_pp contains a pointer to the position in the*/
                       /* string after the token decoded */
  } else {
    in_p=(&in_p[i+1]);
    LenInStr = strlen(in_p);
    i=0;
    while ((in_p[i]!='\"' || in_p[i+1]=='\"') && i<MaxLen && 
           in_p[i]!='\0' && i<LenInStr) {
      if (in_p[i]!='\"') {
        str[i-NQuote]=in_p[i];
        i++;
      } else {
        str[i-NQuote]='\"';
        NQuote++;
        i+=2;
      }
    }
    if (in_p[i]=='\"') {
      *out_pp=(&in_p[i+1]); /*out_pp contains pointer to position in */
                            /* string after the token decoded */
    } else if (i==MaxLen) {
      Error_write(EMAXLENTOK,"gettok",ERR_SEVERE,
                  "Maximum length for token (%d) exceeded. String:\n%s",
                  MaxLen,str);
      return 0;
    } else {
      Error_write(EQUOTE,"gettok",ERR_SEVERE,"Missing closing quote in %s",in_p);
      return 0;
    }
  }

  if (i==MaxLen) {
    Error_write(EMAXLENTOK,"gettok",ERR_SEVERE,
                "Maximum length for token (%d) exceeded in string:\n%s",
                MaxLen,str);
    return 0;
  }

  str[i-NQuote]='\0';

  return 1; 
}

char *StrUtil_substr(char *from, int start, int length) {
  int lenFrom = strlen(from);
  int end;
  char *toString;

  if (start > lenFrom) {
    fprintf(stderr,"Error: substr called with start greater than length of string start=%d str=%s\n",start,from);
    exit(1);
  }
  if (length > 0) {
    end = start + length - 1;
    if (end >= lenFrom) {
      end = lenFrom-1;
    }
  } else if (length == 0) {
    end = start-1;
  } else { // length is negative
    end = lenFrom+length-1;
    if (end < start) {
      fprintf(stderr,"Error: substr called with negative length which resulted in end < start start=%d length=%d str=%s\n",
              start, length, from);
      exit(1);
    }
  }
  return StrUtil_copyNString(&toString, from, start, end-start+1);
}

/******************************************************************************/
/* Routine    :                                                               */
/*             StrUtil_rmspace()                                              */
/* Role       :                                                               */
/*             to remove all white spaces from a string                       */
/* Arguments  :                                                               */
/*             string - the string to remove spaces from                      */
/* Returns    :                                                               */
/*             the length of the stripped string                              */
/******************************************************************************/
int StrUtil_rmspace(char *string) {
  char *chP1 = string;
  char *chP2 = string;
  char ch;
  int inDoubleQuotes = 0;
  int inSingleQuotes = 0;
   
  while ((ch = *chP1) != '\0') {
    if (ch == '\'') inSingleQuotes = (inSingleQuotes) ? 0 : 1;
    if (ch == '"') inDoubleQuotes = (inDoubleQuotes) ? 0 : 1;

    if (inSingleQuotes || inDoubleQuotes || (ch!=' ' && ch!='\t' && ch!='\n')) {
      *chP2=ch;
      chP2++;
    }
    chP1++;
  }
  *chP2='\0';

  return(chP2-string);
}

int StrUtil_truncateAtChar(char *string, char truncCh) {
  char *chP = string;
  char ch;
  int inDoubleQuotes = 0;
  int inSingleQuotes = 0;

  while ((ch = *chP) != '\0') {
    if (ch == '\'') inSingleQuotes = (inSingleQuotes) ? 0 : 1;
    if (ch == '"') inDoubleQuotes = (inDoubleQuotes) ? 0 : 1;

    if (ch==truncCh && !inDoubleQuotes && !inSingleQuotes) {
      break;
    }
    chP++;
  }
  *chP='\0';

  return(chP-string);
}

/******************************************************************************/
/* Routine    :                                                               */
/*             RmQuotes()                                                     */
/* Role       :                                                               */
/*             To remove quotes (if any) from around a string. It is not an   */
/*             error for there to be no quotes.                               */
/* Arguments  :                                                               */
/*             Str - the string to remove quotes from around.                 */
/* Returns    :                                                               */
/*             The length of the new string.                                  */
/* History    :                                                               */
/*             22/03/98 SMJS  Initial implementation                          */
/******************************************************************************/
int StrUtil_rmQuotes(char *str) {
  int len = strlen(str);

  if (str[0]=='\"' || str[0]=='\'') {
    len--;
    memmove(str,&(str[1]),len);
    if (str[len-1] == '\"' || str[len-1] == '\'') {
      str[len-1] = '\0';
      len--;
    }
  }
  return len;
}

int StrUtil_tokenize(char ***tokens, int *ntok, char *string) {
  char token[MAXSTRLEN];
  char *ChP = string;
  int count = 0;

  *tokens = NULL;


  while (StrUtil_gettok(token,&ChP,ChP,MAXSTRLEN)) {
    if (!strlen(token)) break;
    if (!count || !((count+1)%10)) {
      *tokens = (char **)realloc(*tokens,(count+10)*sizeof(char *));
    }
#ifdef DBG
    printf("Adding token %s\n",token);
#endif
    StrUtil_copyString(&((*tokens)[count++]),token,0);
  }
  *ntok = count;
  return 1;
}

int StrUtil_tokenizeByDelim(char ***tokens, int *ntok, char *string, char *delim) {
  char token[MAXSTRLEN];
  char *chP = string;
  char *retChP;
  int count = 0;

  *tokens = NULL;

  while ((retChP = strpbrk(chP,delim))) {
    StrUtil_strncpy(token,chP,retChP-chP);
    if (!count || !((count+1)%10)) {
      *tokens = (char **)realloc(*tokens,(count+10)*sizeof(char *));
    }
    StrUtil_copyString(&((*tokens)[count++]),token,0);
    chP = retChP+1;
    //printf("chP = %c\n",*chP);
  }
// Remainder of string is last token
  if (!count || !((count+1)%10)) {
    *tokens = (char **)realloc(*tokens,(count+10)*sizeof(char *));
  }
  if (strlen(chP)) {
    StrUtil_copyString(&((*tokens)[count++]),chP,0);
  }
  *ntok = count;
  return 1;
}

// Is the entire string a long integer - must be the complete string eg.
//  "1000" - yes
//  "1000 bases" - no
// Returns
//   1 - it is a long int
//   0 - it isn't
// retVal will contain the long value on return if the string is a long int
int StrUtil_isLongInteger(long *retVal, char *str) {
  char *endPtr;
  *retVal = strtol(str, &endPtr, 10);

  if (*endPtr!='\0' || endPtr == str) {
    return 0;
  }
  return 1;
}

/******************************************************************************/
/* Routine    :                                                               */
/*             StrUtil_rmChar()                                               */
/* Role       :                                                               */
/*             To remove all occurences of a character from a string.         */
/* Arguments  :                                                               */
/*             string - the string to remove characters from.                 */
/*             ch     - character to remove                                   */
/* Returns    :                                                               */
/*             the length of the stripped string                              */
/******************************************************************************/
int StrUtil_rmChar(char *string,char ch)
{
   int i=0,j=0;
   int length=0;


   length=strlen(string);
   for (i=0;i<length;i++)
   {
      if (string[i]!=ch)
      {
         string[j]=string[i];
         j++;
      }
   }
   string[j]='\0';
   return(j);
}

int StrUtil_stringCompFunc(const void *a, const void *b) {
  char *s1 = *((char **)a);
  char *s2 = *((char **)b);

  return strcmp(s1,s2);
}

