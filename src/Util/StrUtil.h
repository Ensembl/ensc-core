#ifndef __STRUTIL_H__
#define __STRUTIL_H__

#include <stdio.h>

#define EXTREMELEN 65536
#define MAXSTRLEN 1024

char *StrUtil_copyString(char **to, char *from,int Start);
char *StrUtil_appendString(char *to, char *from);
void StrUtil_reverseString(char *string, int len);
int StrUtil_strupr(char *string);

#ifdef __STRUTIL_MAIN__
 char *emptyString = "";
#else
 extern char *emptyString;
#endif


#endif
