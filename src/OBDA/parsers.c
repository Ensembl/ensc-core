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

#include <string.h>
#include <stdlib.h>

#include "parsers.h"

#include "StrUtil.h"

/* NOTE: tokenParser skips the > */
/* NOTE2: the workStr is fixed length to save having to do a strlen */
/*        and anyway it is highly unlikely that the source of this */
/*        this string read it without using a fixed length buffer. */
/*        MAXSTRLEN is a length I use for fixed length things like */
/*        this. */
/* NOTE3: We need a workStr because strtok mangles any string it */
/*        gets its hands on (deliberately) */

ParserFunc Parser_lookup(char *name) {
  ParserName *curP = parserArray;
 
  while (curP->func != NULL) {
    if (!strcmp(curP->name,name)) {
      return curP->func;
    }
    curP++;
  }

  printf("Error: Didn't find func named %s\n",name);
  return NULL;
}

char *tokenParser(char *header, char *delim, int toknum, char *routine) {
  int i;
  char *token;
  char workStr[MAXSTRLEN];
  char *retToken;


  strcpy(workStr,header);
  token = strtok(workStr+1,delim);
  for (i=0;i<toknum;i++) {
    token = strtok(NULL,delim);
  }
  if (token == NULL) {
    printf("Error parsing ID in %s from string %s\n",routine,header);
  }

  retToken = StrUtil_copyString(&retToken,token,0);

  return retToken;
}

char *prefixParser(char *header, char *delim, char *prefix, char *routine) {
  char *token;
  char workStr[MAXSTRLEN];
  char *retToken;
  int len = strlen(prefix);

  strcpy(workStr,header);
  token = strtok(workStr+1,delim);
  while (token && strncmp(token,prefix,len)) {
    token = strtok(NULL,delim);
  }

  if (token == NULL) {
    printf("Error parsing ID in %s from string %s\n",routine,header);
  }

  retToken = StrUtil_copyString(&retToken,token,len);

  return retToken;
}

char *multiTokenParser(char *header, char *delim, int toknum,char **terminators,
                       int nterminator, char *routine, char **retPos) {
  int i;
  char *token;
  char workStr[MAXSTRLEN];
  char *retToken;

/*
  printf("Got header %s\n",header);
*/
  strcpy(workStr,header);
  token = strtok(workStr,delim);
  for (i=0;i<toknum;i++) {
    token = strtok(NULL,delim);
  }


  if (token == NULL ) {
    return NULL;
  }

  for (i=0;i<nterminator;i++) {
    if (!strcmp(token,terminators[i])) {
      return NULL;
    }
  }

  *retPos = header + (token - workStr) + strlen(token); 
   
  retToken = StrUtil_copyString(&retToken,token,0);

  return retToken;
}

char *unigeneGenbankParser(char *header) {
  return prefixParser(header," ","/gb=","unigeneGenbankParser");
}

char *unigeneParser(char *header) {
  return prefixParser(header," ","/ug=","unigeneParser");
}

char *dbParser(char *header) {
  return tokenParser(header,"|",3,"dbParser");
}

char *dbNoVersionParser(char *header) {
  char *id = tokenParser(header,"|",3,"dbNoVersionParser");
  char *ChP;
  if ((ChP = strchr(id,'.')) != NULL) {
    *ChP = '\0';
  }
  return id;
}

char *rikenParser(char *header) {
  return tokenParser(header,"|",2,"rikenParser");
}

/* This is unused */
char *dbIDParser(char *header) {
  return tokenParser(header,"|",1,"dbIDParser");
}

char *emblParser(char *header) {
  return tokenParser(header," \t",1,"emblParser");
}

char *locuslinkParser(char *header) {
  return tokenParser(header," ",0,"locuslinkParser");
}

char *singleWordParser(char *header) {
  return tokenParser(header," \n",0,"singleWordParser");
}

/* NOTE: Not sure whether . should be included */
char *refseqIDParser(char *header) {
  return tokenParser(header,"|",1,"refseqIDParser");
}

/* NOTE: Not sure whether . should be included */
char *refseqParser(char *header) {
  return tokenParser(header,"|",3,"refseqParser");
  return NULL;
}

char *swallIDParser(char *header) {
  return tokenParser(header," \t",0,"swallIDParser");
}

char *uniprotIDParser(char *header) {
  return tokenParser(header,".\t \n",0,"uniprotIDParser");
}

char *swallParser(char *header) {
  return tokenParser(header," \t",1,"swallParser");
}

char *swissIDParser(char *header) {
  return tokenParser(header," \t",0,"swissIDParser");
}

char *swissParser(char *header) {
  return tokenParser(header," \t",1,"swissParser");
}

/* NOTE: This parser uses : as the delimiter */
char *swirIDParser(char *header) {
  return tokenParser(header,": ",1,"swirIDParser");
}

char *swirParser(char *header) {
  return tokenParser(header," \t",1,"swirParser");
}

char *traceParser(char *header) {
  return tokenParser(header," \t",0,"traceParser");
}

char *tremblParser(char *header) {
  return tokenParser(header," \t",0,"tremblParser");
}
char *rodParser(char *header) {
  return tokenParser(header," \t\n",1,"rodParser");
}
char *rodIDParser(char *header) {
  return tokenParser(header," \t",1,"rodIDParser");
}
char *wormParser(char *header) {
  return tokenParser(header," \t\n;",1,"wormParser");
}
char *wormIDParser(char *header) {
  return tokenParser(header," \t",1,"wormIDParser");
}

char *refseqProt2RNAParser(char *header) {
  char *brastr = tokenParser(header," ",1,"refseqProt2RNAParser");
  int len = strlen(brastr);
  int i; 
  if (brastr[0] == '(' && brastr[len-1] == ')') {
    for (i=1;i<len-1;i++) {
      brastr[i-1] = brastr[i];
    }
    brastr[i-1] = '\0';
  } else {
    printf("Refseq header doesn't match expected format %s\n",header);
  }
  return brastr;
}

char *swallMultiParser(char *header, char **retPos) {
  static char *terminators[] = {"Desc:","-"};
  char *chP;
  char *id;

  if (*retPos == NULL) {
/* first call */
    id = multiTokenParser(header," \t\n",2,terminators,2,
                          "swallMultiParser",retPos);
  } else {
/* subsequent calls */
    id = multiTokenParser(header," \t\n",0,terminators,2,
                          "swallMultiParser",retPos);
  }
  if (id && (chP = strchr(id,'.'))) {
    *chP = '\0';
  }
  return id;
}

char *BTMultiParser(char *header, char **retPos) {
  static char *terminators[] = {};
  char *id;

  if (*retPos == NULL) {
/* first call */
    id = multiTokenParser(header," \t\n",1,terminators,0,
                          "BTMultiParser",retPos);
  } else {
/* subsequent calls */
    id = multiTokenParser(header," \t\n",0,terminators,0,
                          "BTMultiParser",retPos);
  }
  return id;
}
