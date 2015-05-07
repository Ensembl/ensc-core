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

/****************************************************************\
*                                                                *
*  Library of parsing routines for FASTA format sequence data    *
*                                                                *
*  Steve Searle         mailto:searle@sanger.ac.uk               *
*  Copyright (C) 2002.  All Rights Reserved.                     *
*                                                                *
*  This source code is distributed under the terms of the        *
*  GNU Lesser General Public License. See the file COPYING       *
*  or http://www.fsf.org/copyleft/lesser.html for details        *
*                                                                *
*  If you use this code, please keep this notice intact.         *
*                                                                *
\****************************************************************/

#ifndef INCLUDED_PARSERS_H
#define INCLUDED_PARSERS_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <stdio.h>

char *unigeneGenbankParser(char *header);
char *unigeneParser(char *header);
char *dbParser(char *header);
char *dbNoVersionParser(char *header);
char *rikenParser(char *header);
char *dbIDParser(char *header);
char *emblParser(char *header);
char *locuslinkParser(char *header);
char *singleWordParser(char *header);
char *refseqIDParser(char *header);
char *refseqParser(char *header);
char *swallIDParser(char *header);
char *swallParser(char *header);
char *swissIDParser(char *header);
char *swissParser(char *header);
char *swirIDParser(char *header);
char *swirParser(char *header);
char *traceParser(char *header);
char *tremblParser(char *header);
char *rodParser(char *header);
char *rodIDParser(char *header);
char *uniprotIDParser(char *header);
char *wormParser(char *header);
char *wormIDParser(char *header);
char *refseqProt2RNAParser(char *header);
char *swallMultiParser(char *header, char **retPos);
char *BTMultiParser(char *header, char **retPos);

typedef char * (*ParserFunc)(char *);

ParserFunc Parser_lookup(char *name);

typedef struct ParserNameStruct {
  ParserFunc func;
  char *name;
} ParserName;

static ParserName parserArray[] = {
  { unigeneGenbankParser, "unigeneGenbankParser" },
  { unigeneParser, "unigeneParser" },
  { dbParser, "dbParser" },
  { dbNoVersionParser, "dbNoVersionParser" },
  { rikenParser, "rikenParser" },
  { dbIDParser, "dbIDParser" },
  { emblParser, "emblParser" },
  { locuslinkParser, "locuslinkParser" },
  { singleWordParser, "singleWordParser" },
  { refseqIDParser, "refseqIDParser" },
  { refseqParser, "refseqParser" },
  { swallIDParser, "swallIDParser" },
  { swallParser, "swallParser" },
  { swissIDParser, "swissIDParser" },
  { swissParser, "swissParser" },
  { swirIDParser, "swirIDParser" },
  { swirParser, "swirParser" },
  { traceParser, "traceParser" },
  { tremblParser, "tremblParser" },
  { rodParser, "rodParser" },
  { rodIDParser, "rodIDParser" },
  { wormParser, "wormParser" },
  { wormIDParser, "wormIDParser" },
  { refseqProt2RNAParser, "refseqProt2RNAParser" },
  { swallMultiParser, "swallMultiParser" },
  { BTMultiParser, "BTMultiParser" },
  { uniprotIDParser, "uniprotIDParser" },
  { NULL, NULL }
};


#define MAXSTRLEN 1024

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_PARSERS_H */
