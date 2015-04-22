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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "CHash.h"
#include "Error.h"
#include "StrUtil.h"
#include "Stream.h"

/******************************************************************************/
/* Routine    :                                                               */
/*             CHash_addAllocedStr()                                          */
/* Role       :                                                               */
/*             To add an extra string to a CHASHTABLE by reallocing the array */
/* Arguments  :                                                               */
/*             Table    - the created hash table (RETURNED CHANGED)           */
/*             String   - the string to add                                   */
/* Returns    :                                                               */
/*             1 - success                                                    */
/*             0 - failure                                                    */
/* History    :                                                               */
/*             31/08/98 SMJS  Initial Implementation                          */
/******************************************************************************/
int CHash_addAllocedStr(CHASHTABLE *Table,char *String) {

  if (String==NULL || Table==NULL) {
    Error_write(ENULLPOINTER,"CHash_addAllocedStr",ERR_SEVERE,"String or Table was NULL");
    return 0;
  }

  if (!Table->NElement || !(Table->NElement%CHASH_ALLOCSIZE)) {
    if ((Table->Strings = (char **)realloc(Table->Strings,
                                           (Table->NElement+CHASH_ALLOCSIZE)*
                                           sizeof(char *)))==NULL) {
      Error_write(EREALLERR,"CHash_addAllocedStr",ERR_SEVERE,"Table->Strings");
      return 0;
    }
  }

#ifdef DBG
  Stream_fprintf(DBGStream,"Adding %s\n",String);
#endif

  Table->Strings[Table->NElement] = String;

  if (!CHash_insert(String,Table,Table->NElement)) {
    Error_trace("CHash_addAllocedStr",String);
    return 0;
  }

  Table->NElement++;
   
  return 1;
}

/******************************************************************************/
/* Routine    :                                                               */
/*             CHash_addStr()                                                 */
/* Role       :                                                               */
/*             To add an extra string to a CHASHTABLE by reallocing the array */
/* Arguments  :                                                               */
/*             Table    - the created hash table (RETURNED CHANGED)           */
/*             String   - the string to add                                   */
/* Returns    :                                                               */
/*             1 - success                                                    */
/*             0 - failure                                                    */
/* History    :                                                               */
/*             31/08/98 SMJS  Initial Implementation                          */
/******************************************************************************/
int CHash_addStr(CHASHTABLE *Table,char *String) {
  if (String==NULL) {
    Error_write(ENULLPOINTER,"CHash_addStr",ERR_SEVERE, "String was NULL");
    return 0;
  }

  if (!Table->NElement || !(Table->NElement%CHASH_ALLOCSIZE)) {
    if ((Table->Strings=(char **)realloc(Table->Strings,
                                         (Table->NElement+CHASH_ALLOCSIZE)*
                                          sizeof(char *)))==NULL) {
      Error_write(EREALLERR,"CHash_addStr",ERR_SEVERE,"Table->Strings");
      return 0;
    }
  }

#ifdef DBG
  Stream_fprintf(DBGStream,"Adding %s\n",String);
#endif

  if (!StrUtil_copyString(&(Table->Strings[Table->NElement]),String,0)) {
    Error_trace("CHash_addStr",String);
    return 0;
  }

  if (!CHash_insert(String,Table,Table->NElement)) {
    Error_trace("CHash_addStr",String);
    return 0;
  }

  Table->NElement++;
   
  return 1;
}

/******************************************************************************/
/* Routine    :                                                               */
/*             CHash_addToArray()                                             */
/* Role       :                                                               */
/*             Add a string to a CHASHARRAY at a specified position           */
/* Arguments  :                                                               */
/*             String   - the array containing the unalphebetized strings     */
/*             Table    - the number of elements in Array                     */
/* Returns    :                                                               */
/*             1 - success                                                    */
/*             0 - failure                                                    */
/* History    :                                                               */
/*             28/02/97 SMJS  Initial Implementation                          */
/******************************************************************************/
int CHash_addToArray(CHASHTABLE *Table,CHASHARRAY *LetArray,int NArray,
                    char *String,int StrInd, int Position) {
  int i = 0;
  if (LetArray==NULL) {
    Error_write(ENULLPOINTER,"CHash_addStr",ERR_SEVERE, "LetArray was NULL");
    return 0;
  }
   
  for (i=NArray;i>Position;i--) {
    LetArray[i].Index=LetArray[i-1].Index;
  }

  LetArray[Position].Index=StrInd;
   
  return 1;
}

/******************************************************************************/
/* Routine    :                                                               */
/*             CHash_binSearch()                                              */
/* Role       :                                                               */
/*             To search for a string in a CHASHARRAY                         */
/* Arguments  :                                                               */
/*             String   - the array containing the unalphebetized strings     */
/*             Table    - the number of elements in Array                     */
/* Returns    :                                                               */
/*             The index of the matching entry in the CHASHARRAY OR           */
/*             -1 if not found.                                               */
/* History    :                                                               */
/*             28/02/97 SMJS  Initial Implementation                          */
/******************************************************************************/
int CHash_binSearch(CHASHTABLE *Table,CHASHARRAY *Array,char *String,int *Ind, 
                    int low, int high) {
  int middle = 0;
  if (Table==NULL || Array==NULL || String==NULL || Ind==NULL) {
    Error_write(ENULLPOINTER,"CHash_addStr",ERR_SEVERE, "Table, Array, String or Ind was NULL");
    return 0;
  }

  if (low>high) {
    *Ind = low;
    return -1;
  } else {
    middle=(high+low)/2;
    if (!strcmp(String,Table->Strings[Array[middle].Index])) {
      *Ind = middle;
      return middle;
    } else if (strcmp(String,Table->Strings[Array[middle].Index])<0) {
      return CHash_binSearch(Table,Array,String,Ind,low,middle-1);
    } else {
      return CHash_binSearch(Table,Array,String,Ind,middle+1,high);
    }
  }
}

/******************************************************************************/
/* Routine    :                                                               */
/*             CHash_dump()                                                   */
/* Role       :                                                               */
/*             To dump a CHASHTABLE for debugging                             */
/* Arguments  :                                                               */
/*             Table    - the CHASHTABLE to dump                              */
/* Returns    :                                                               */
/*             None                                                           */
/* History    :                                                               */
/*             02/03/97 SMJS  Initial Implementation                          */
/******************************************************************************/
void CHash_dump(CHASHTABLE *Table) {
  int i;
  int j;

  if (Table==NULL) {
    Error_write(ENULLPOINTER,"CHash_addStr",ERR_SEVERE, "Table was NULL");
    return;
  }

  Stream_fprintf(DBGStream,"DBG>> Dumping CHash Table\n      ===================\n");
  Stream_fprintf(DBGStream,"DBG>> NElement = %d\n",Table->NElement);
  for (i=0;i<Table->NElement;i++) {
    if (Table->Strings[i]!=NULL) {
      Stream_fprintf(DBGStream,"DBG>> Strings[%d] = %s\n",i,
      Table->Strings[i]);
    } else {
      Stream_fprintf(DBGStream,"DBG>> Strings[%d] DELETED\n",i);
    }
  }

  Stream_fprintf(DBGStream,"DBG>> \n");

  for (i=0;i<Table->NLetter;i++) {
    if (isprint(i) && Table->HashFunc == CHash_getLetInd) {
      Stream_fprintf(DBGStream,"DBG>> Letter %d (%c) (%d entries)\n",
              i,i/*+'A'*/,Table->LetCounts[i]);
    } else {
      Stream_fprintf(DBGStream,"DBG>> Letter %d (%d entries)\n",
              i,Table->LetCounts[i]);
    }
    for (j=0;j<Table->LetCounts[i];j++) {
      Stream_fprintf(DBGStream,"DBG>> Letter[%d][%d].Index = %d (%s)\n",
              i,j,Table->Letter[i][j].Index,
      Table->Strings[Table->Letter[i][j].Index]);
    }
    Stream_fprintf(DBGStream,"DBG>> \n");
  }

  return;
}

/******************************************************************************/
/* Routine    :                                                               */
/*             CHash_find()                                                   */
/* Role       :                                                               */
/*             To find the index of a string in a CHASHTABLE                  */
/* Arguments  :                                                               */
/*             String   - the string to find                                  */
/*             Table    - the number of elements in Array                     */
/*             StrInd   - the index of the string (RETURNED)                  */
/* Returns    :                                                               */
/*             1 - success                                                    */
/*             0 - failure                                                    */
/* History    :                                                               */
/*             25/02/97 SMJS  Initial Implementation                          */
/* Note       :                                                               */
/*             No error is given if the string is not in the array            */
/******************************************************************************/
int CHash_find(char *String,CHASHTABLE *Table,int *StrInd) {
  int LetInd;
  int ArrayLoc = 0;
   
  if (Table==NULL) {
    Error_write(ENULLPOINTER,"CHash_addStr",ERR_SEVERE, "Table was NULL");
    return 0;
  }

  if (!Table->HashFunc(String,&LetInd)) {
    Error_trace("CHash_find",String);
    return 0;
  }
   
  *StrInd = -1;
  if (Table->LetCounts[LetInd]) {
    if (CHash_binSearch(Table,Table->Letter[LetInd],String,&ArrayLoc,0,
                       Table->LetCounts[LetInd]-1) == -1) {
      return 0;
    }
    *StrInd = Table->Letter[LetInd][ArrayLoc].Index;
    return 1;
  }

  return 0;
}

/******************************************************************************/
/* Routine    :                                                               */
/*             CHash_firstFour()                                              */
/* Role       :                                                               */
/*             To get letter index of String using first four letters         */
/* Arguments  :                                                               */
/*             String   - the string to get index for                         */
/*             LetInd   - index (RETURNED)                                    */
/* Returns    :                                                               */
/*             1 - success                                                    */
/*             0 - failure                                                    */
/* History    :                                                               */
/*             07/03/00 SMJS  Initial Implementation                          */
/******************************************************************************/
int CHash_firstFour(char *String,int *LetInd) {
  int   i;
  char *ChP;
   
  if (String==NULL || LetInd==NULL) {
    Error_write(ENULLPOINTER,"CHash_addStr",ERR_SEVERE, "String or LetInd was NULL");
    return 0;
  }

  *LetInd=0;
  for (i=0,ChP=String;i<4 && *ChP!='\0';i++,ChP++) {
    *LetInd += *ChP;
  }
  return 1;
}

/******************************************************************************/
/* Routine    :                                                               */
/*             CHash_fourLets()                                               */
/* Role       :                                                               */
/*             To get letter index of String                                  */
/* Arguments  :                                                               */
/*             String   - the string to get index for                         */
/*             LetInd   - index (RETURNED)                                    */
/* Returns    :                                                               */
/*             1 - success                                                    */
/*             0 - failure                                                    */
/* History    :                                                               */
/*             04/09/98 SMJS  Initial Implementation                          */
/******************************************************************************/
int CHash_fourLets(char *String,int *LetInd) {
  int   i;
  char *ChP;

  if (String==NULL || LetInd==NULL) {
    Error_write(ENULLPOINTER,"CHash_addStr",ERR_SEVERE, "String or LetInd was NULL");
    return 0;
  }

  *LetInd=0;
  for (i=0,ChP=String;i<10 && *ChP!='\0';i++,ChP++) {
    if (i==0 || i==3 || i==6 || i==9) {
      *LetInd += *ChP;
    }
  }
  return 1;
}

/******************************************************************************/
/* Routine    :                                                               */
/*             CHash_free()                                                   */
/* Role       :                                                               */
/*             To free a CHASHTABLE                                           */
/* Arguments  :                                                               */
/*             Table    - the table                                           */
/* Returns    :                                                               */
/*             None                                                           */
/* History    :                                                               */
/*             02/03/97 SMJS Initial Implementation                           */
/******************************************************************************/
void CHash_free(CHASHTABLE *Table) {
  int i;

  if (Table==NULL) {
    Error_write(ENULLPOINTER,"CHash_addStr",ERR_SEVERE, "Table was NULL");
    return;
  }

  for (i=0;i<Table->NElement;i++) {
    if (Table->Strings[i]) {
      free(Table->Strings[i]);
    }
  }
  for (i=0;i<Table->NLetter;i++) {
    if (Table->Letter[i]) {
      free(Table->Letter[i]);
    }
  }
  free(Table->Letter);
  free(Table->LetCounts);
  free(Table);
  return;
}

/******************************************************************************/
/* Routine    :                                                               */
/*             CHash_getLetInd()                                              */
/* Role       :                                                               */
/*             To get letter index of String                                  */
/* Arguments  :                                                               */
/*             String   - the string to get index for                         */
/*             LetInd   - index (RETURNED)                                    */
/* Returns    :                                                               */
/*             1 - success                                                    */
/*             0 - failure                                                    */
/* History    :                                                               */
/*             28/02/97 SMJS  Initial Implementation                          */
/******************************************************************************/
int CHash_getLetInd(char *String,int *LetInd) {

  if (String==NULL) {
    Error_write(ENULLPOINTER,"CHash_addStr",ERR_SEVERE, "String was NULL");
    return 0;
  }

  *LetInd = String[0];

  return 1;
}

/******************************************************************************/
/* Routine    :                                                               */
/*             CHash_getSorted()                                              */
/* Role       :                                                               */
/*             Get an array of string pointers to the entries in a CHASH table*/
/*             in sorted order.                                               */
/* Arguments  :                                                               */
/*             Table    - the CHASH table                                     */
/*             Array    - the array containing the sorted string pointers     */
/*             NArray   - the number of elements in Array                     */
/* Returns    :                                                               */
/*             1 - success                                                    */
/*             0 - failure                                                    */
/* History    :                                                               */
/*             26/08/98 SMJS  Initial Implementation                          */
/******************************************************************************/
int CHash_getSorted(CHASHTABLE *Table,char ***Array,int *NArray) {
  int i;
  int j;
  int NAdded = 0;
   
  if (Table==NULL) {
    Error_write(ENULLPOINTER,"CHash_addStr",ERR_SEVERE, "Table was NULL");
    return 0;
  }

  if (!Table->NElement) {
    /* Error */
    return 0;
  }
  if ((*Array = (char **)calloc(Table->NElement,sizeof(char *)))==NULL) {
    Error_write(ECALLERR,"CHash_getSorted",ERR_SEVERE,"*Array");
    return 0;
  }
   
  for (i=0;i<Table->NLetter;i++) {
    for (j=0;j<Table->LetCounts[i];j++) {
      (*Array)[NAdded++] = Table->Strings[Table->Letter[i][j].Index];
    }
  }

/* NAdded isn't necessarily the same as NElement because of deletions from a */
/* CHASHTABLE */
  *NArray = NAdded;

  return 1;
}

/******************************************************************************/
/* Routine    :                                                               */
/*             CHash_init()                                                   */
/* Role       :                                                               */
/*             To initialise a hash table. This involves allocating the table,*/
/*             setting the hashing function and allocating the array          */
/* Arguments  :                                                               */
/*             Table    - the created hash table (RETURNED ALLOCATED)         */
/*             Func     - the hash function to use in this hash table.        */
/*             NLetter  - the number of elements in Func                      */
/* Returns    :                                                               */
/*             1 - success                                                    */
/*             0 - failure                                                    */
/* History    :                                                               */
/*             04/09/98 SMJS  Initial Implementation                          */
/******************************************************************************/
int CHash_init(CHASHTABLE **Table,int Func(),int NLetter) {

  if (!NLetter) {
    Error_write(ECHASH,"CHash_init",ERR_SEVERE,"Number of buckets is Zero");
    return 0;
  }
  if ((*Table = (CHASHTABLE *)calloc(1,sizeof(CHASHTABLE)))==NULL) {
    Error_write(ECALLERR,"CHash_init",ERR_SEVERE,"*Table");
    return 0;
  }

  (*Table)->HashFunc = Func;
  (*Table)->NLetter  = NLetter;

  if (((*Table)->Letter=(CHASHARRAY **)calloc(NLetter,sizeof(CHASHARRAY *)))==NULL) {
    Error_write(ECALLERR,"CHash_init",ERR_SEVERE,"*Table->Letter");
    return 0;
  }

  if (((*Table)->LetCounts=(int *)calloc(NLetter,sizeof(int)))==NULL) {
    Error_write(ECALLERR,"CHash_init",ERR_SEVERE,"*Table->LetCounts");
    return 0;
  }
  return 1;
}

/******************************************************************************/
/* Routine    :                                                               */
/*             CHash_insert()                                                 */
/* Role       :                                                               */
/*             To insert a string into a CHASHTABLE                           */
/* Arguments  :                                                               */
/*             String   - the array containing the unalphebetized strings     */
/*             Table    - the number of elements in Array                     */
/* Returns    :                                                               */
/*             1 - success                                                    */
/*             0 - failure                                                    */
/* History    :                                                               */
/*             25/02/97 SMJS  Initial Implementation                          */
/******************************************************************************/
int CHash_insert(char *String,CHASHTABLE *Table,int StrInd) {
  int LetInd;
  int ArrayLoc=0;

  if (Table==NULL) {
    Error_write(ENULLPOINTER,"CHash_addStr",ERR_SEVERE, "Table was NULL");
    return 0;
  }

  /* First find which array to add to */
  if (!Table->HashFunc(String,&LetInd)) {
    Error_trace("CHash_insert",String);
    return 0;
  }

  /* Now check if we need to (re)allocate it */
  if (!Table->LetCounts[LetInd]) {
    if ((Table->Letter[LetInd]=
        (CHASHARRAY *)calloc(10,sizeof(CHASHARRAY)))==NULL) {
      Error_write(ECALLERR,"CHash_insert",ERR_SEVERE,"Table->Letter");
      return 0;
    }
  } else if (!(Table->LetCounts[LetInd]%10)) {
    if ((Table->Letter[LetInd]=
        (CHASHARRAY *)realloc(Table->Letter[LetInd],
         (Table->LetCounts[LetInd]+10)*sizeof(CHASHARRAY)))==NULL) {
      Error_write(EREALLERR,"CHash_insert",ERR_SEVERE,"Table->Letter");
      return 0;
    }
  }

/* Find the location in the array to add the string */
  if (Table->LetCounts[LetInd]) {
    if (CHash_binSearch(Table,Table->Letter[LetInd],String,&ArrayLoc,0,
                       Table->LetCounts[LetInd]-1)!= -1) {
      Error_trace("CHash_insert",String);
      return 0;
    }
  }

/* Now add the string in its sorted location */
  if (!CHash_addToArray(Table,Table->Letter[LetInd],Table->LetCounts[LetInd],
                       String,StrInd,ArrayLoc)) {
    Error_trace("CHash_insert",String);
    return 0;
  }
  (Table->LetCounts[LetInd])++;
  return 1;
}

/******************************************************************************/
/* Routine    :                                                               */
/*             CHash_removeStr()                                              */
/* Role       :                                                               */
/*             To remove a string from a CHASHTABLE                           */
/* Arguments  :                                                               */
/*             Table    - the hash table (RETURNED CHANGED)                   */
/*             String   - the string to remove                                */
/*             StrInd   - the index of String in Table->Strings               */
/*             LetInd   - the index of Letter for string in Table->Letters    */
/* Returns    :                                                               */
/*             1 - success                                                    */
/*             0 - failure                                                    */
/* History    :                                                               */
/*             01/09/98 SMJS  Initial Implementation                          */
/******************************************************************************/
int CHash_removeStr(CHASHTABLE *Table,char *String,int StrInd,int LetInd) {
  int i;
  int LetPos;

  if (Table==NULL || String==NULL) {
    Error_write(ENULLPOINTER,"CHash_addStr",ERR_SEVERE, "Table or String was NULL");
    return 0;
  }
   
#ifdef DBG
  Stream_fprintf(DBGStream,"Removing %s\n",String);
#endif
   
/* Find index in Letter[LetInd] */
  if (Table->Letter[LetInd] == NULL) {
    Error_write(ECHASH,"CHash_removeStr",ERR_SEVERE,"Table->Letter[%d] is NULL",
             LetInd);
    return 0;
  }

  for (LetPos=0;LetPos<Table->LetCounts[LetInd];LetPos++) {
    if (Table->Letter[LetInd][LetPos].Index == StrInd)
      break;
  }

  if (LetPos==Table->LetCounts[LetInd]) {
    Error_write(ECHASH,"CHash_removeStr",ERR_SEVERE,
             "Didn't find Index %d in Table->Letter[%d]",StrInd,LetInd);
    return 0;
  }

  for (i=LetPos;i<Table->LetCounts[LetInd]-1;i++) {
    Table->Letter[LetInd][i].Index=Table->Letter[LetInd][i+1].Index;
  }

  (Table->LetCounts[LetInd])--;

  if (!Table->LetCounts[LetInd])
    free(Table->Letter[LetInd]);
   
  free(Table->Strings[StrInd]);
  Table->Strings[StrInd] = NULL;
   
  return 1;
}

/******************************************************************************/
/* Routine    :                                                               */
/*             CHash_setup()                                                  */
/* Role       :                                                               */
/*             To setup a simple alphabetic hash table                        */
/* Arguments  :                                                               */
/*             Array    - the array containing the unalphebetized strings     */
/*             NArray   - the number of elements in Array                     */
/*             Table    - the created hash table (RETURNED ALLOCATED)         */
/* Returns    :                                                               */
/*             1 - success                                                    */
/*             0 - failure                                                    */
/* History    :                                                               */
/*             ??/??/97 SMJS  Initial Implementation                          */
/*             30/08/98 SMJS  Added ing check                                 */
/*             31/08/98 SMJS  Now allows allocation without adding strings    */
/******************************************************************************/
int CHash_setup(char **Array,int NArray,CHASHTABLE **Table,int Func(),int NLet) {
  int i;

  if (!CHash_init(Table,Func,NLet)) {
    Error_trace("CHash_setup",NULL);
    return 0;
  }

  if (NArray!=0) {
    if (((*Table)->Strings=(char **)calloc(NArray,sizeof(char *)))==NULL) {
      Error_write(ECALLERR,"CHash_setup",ERR_SEVERE,"*Table->Strings");
      return 0;
    }
    for (i=0;i<NArray;i++) {
#ifdef DBG
      Stream_fprintf(DBGStream,"Adding %s\n",Array[i]);
#endif
      if (!StrUtil_copyString(&((*Table)->Strings[i]),Array[i],0)) {
        Error_trace("CHash_setup",Array[i]);
        return 0;
      }
      if (!CHash_insert(Array[i],*Table,i)) {
        Error_trace("CHash_setup",Array[i]);
        return 0;
      }
      (*Table)->NElement++;
    }
  }
  return 1;
}
