#ifndef __STRUTIL_H__
#define __STRUTIL_H__

#include <stdio.h>

#define EXTREMELEN 65536
#define MAXSTRLEN 1024

char *StrUtil_copyString(char **to, char *from,int Start);
char *StrUtil_appendString(char *to, char *from);
void StrUtil_reverseString(char *string, int len);
int StrUtil_strupr(char *string);
int StrUtil_strlwr(char *string);
int StrUtil_rmTrailSpace(char *str);
int StrUtil_rmLeadSpace(char *str);
int StrUtil_rmEndSpaces(char *str);
char *StrUtil_strncpy(char *s1,char *s2,int len);
char *StrUtil_fastCopyString(char **to, char *from,int length);
char *StrUtil_copyNString(char **to, char *from,int Start,int length);
char *StrUtil_copyNTString(char **to, char *from,int Start,int length);
char *StrUtil_copyTString(char **to, char *from,int Start);
int StrUtil_stripnewline(char *str);
int StrUtil_gettok(char *str,char **out_pp,char *in_p,int MaxLen);












#ifdef __STRUTIL_MAIN__
 char *emptyString = "";
#else
 extern char *emptyString;
#endif


#endif
