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
