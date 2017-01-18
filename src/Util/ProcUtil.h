/*
 * Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
 * Copyright [2016-2017] EMBL-European Bioinformatics Institute
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

#ifndef __PROCUTIL_H__
#define __PROCUTIL_H__

#ifdef linux
#define UNW_LOCAL_ONLY
//#include "libunwind.h"
#endif

int ProcUtil_timeInfo(char *routine);
void ProcUtil_mallInfo(void);

#if 0
void ProcUtil_showBacktrace(char *prog);
#ifdef linux
int ProcUtil_getFileAndLine(char *prog, unw_word_t addr, char *file, size_t flen, int *line);
#endif
#endif

#endif
