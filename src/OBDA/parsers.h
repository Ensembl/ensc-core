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
char *wormParser(char *header);
char *wormIDParser(char *header);
char *refseqProt2RNAParser(char *header);
char *swallMultiParser(char *header, char **retPos);

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
  { NULL, NULL }
};


#define MAXSTRLEN 1024

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INCLUDED_PARSERS_H */
