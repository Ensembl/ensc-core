#ifndef __FILEROUT_H__
#define __FILEROUT_H__

void        FileUtil_close(FILE *fp,char *fname);
int         FileUtil_countLines(FILE *fp,char *string);
int         FileUtil_countToLine(FILE *fp,int *NLine,char *string);
int         FileUtil_nextPrefixedLine(FILE *fp,char *string,char *line);
int         FileUtil_exists(char *FName);
int         FileUtil_getRealPath(char *FNIn, char *Path);
int         FileUtil_getStrippedLine(char *Line,int MaxLen, FILE *Fp);
FILE       *FileUtil_open(char *fname,char *stat,char *routine);
int         FileUtil_readToChar(FILE *Fp,char *String,int StrLen,char EndChar);
int         FileUtil_stripExt(char *fname,char *result);
int         FileUtil_stripPath(char *fname,char *result);

#endif /* __FILEROUT_H__ */
