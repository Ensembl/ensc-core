#include <ctype.h>
#include <string.h>

#define __STRUTIL_MAIN__
#include "StrUtil.h"
#undef __STRUTIL_MAIN__
#include "Error.h"

char *StrUtil_appendString(char *to, char *from) {
  int lenTo = strlen(to);
  int lenFrom = strlen(from);

  if ((to = (char *)realloc(to,lenTo+lenFrom+1))== NULL) {
    fprintf(stderr,"ERROR: Failed reallocating string\n");
    exit(1);
  }

  strcat(to,from);
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
    len++;
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
      Error_write(EMALLERR,"CopyString",ERR_SEVERE,"to");
      return NULL;
    }
    (*to)[0]='\0';
  } else {
    Len-=Start;
    cp = &from[Start];
    if ((*to = (char *)calloc(Len+2,sizeof(char)))==NULL) {
      Error_write(EMALLERR,"CopyString",ERR_SEVERE,"to");
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
  int lenTo;
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
  int len=0;


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

