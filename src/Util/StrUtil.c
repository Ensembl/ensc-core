#include "StrUtil.h"
#include <ctype.h>

char *StrUtil_CopyString(char *str) {
  char *copy;
  int len = strlen(str);

  if ((copy = (char *)malloc(len+1))== NULL) {
    fprintf(stderr,"ERROR: Failed allocating string\n");
    return NULL;
  }

  strncpy(copy,str,len);
  copy[len] = '\0';

  return copy;
}

void StrUtil_ReverseString(char *string, int len) {
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

int StrUtil_strupr(char *string) {
  int len=0;;
  char *chP = string;

  while (*chP != '\0') {
    *chP = toupper(*chP);
    len++;
    chP++;
  }
  return len;
}
