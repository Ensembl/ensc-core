/*
 * Copyright [1999-2017] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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


#include "EcoString.h"
#include "Error.h"
#include "CHash.h"
#include "StrUtil.h"
/******************************************************************************/
/* Routine    :                                                               */
/*             EcoString_addStr()                                             */
/* Role       :                                                               */
/*             Add a string to an ECOSTRTABLE.                                */
/* Arguments  :                                                               */
/*             EcoSTabP - Pointer to ECOSTRTABLE to add to                    */
/*             String   - String to add                                       */
/*             StrInd   - Index of string in ECOSTRTABLE (RETURNED)           */
/* Returns    :                                                               */
/*             1 - success                                                    */
/*             0 - failure                                                    */
/* History    :                                                               */
/*             31/08/98 SMJS  Initial Implementation                          */
/******************************************************************************/
int EcoString_addStr(ECOSTRTABLE *EcoSTabP, char *String, int *StrInd) {
  int NElement;

  if (!CHash_addAllocedStr(EcoSTabP->CHashTab,String)) {
    Error_trace("EcoString_addStr",NULL);
    return 0;
  }

  NElement = EcoSTabP->CHashTab->NElement;

  if (NElement==1 || !(NElement%10)) {
    if ((EcoSTabP->UseCount = (int *)realloc(EcoSTabP->UseCount,
                                             (NElement+10)*sizeof(int)))==NULL) {
      Error_write(EREALLERR,"EcoString_addStr",ERR_SEVERE, "EcoSTabP->UseCount");
      return 0;
    }
  
    if ((EcoSTabP->CheckSums = (int *)realloc(EcoSTabP->CheckSums,
                                              (NElement+10)*sizeof(int)))==NULL) {
      Error_write(EREALLERR,"EcoString_addStr",ERR_SEVERE,"EcoSTabP->UseCount");
      return 0;
    }
  }

  if (!EcoSTabP->UseCount || !EcoSTabP->CheckSums) {
    Error_write(EREALLERR,"EcoString_addStr",ERR_SEVERE,"EcoSTabP->UseCount");
    return 0;
  }

  EcoSTabP->UseCount[NElement-1] = 1;
  EcoSTabP->CheckSums[NElement-1] = EcoString_genCheckSum(String);
  *StrInd = NElement-1;

/* Return success */
  return 1;
}

/******************************************************************************/
/* Routine    :                                                               */
/*             EcoString_changeStr()                                          */
/* Role       :                                                               */
/*             Get a pointer to a (possibly new) EcoString                    */
/* Arguments  :                                                               */
/*             EcoSTabP - Pointer to ECOSTRTABLE to add to                    */
/*             To       - Pointer to ECOSTRING (RETURNED)                     */
/*             From     - String to copy from                                 */
/*             StartPos - the first position in From to copy                  */
/* Returns    :                                                               */
/*             1 - success                                                    */
/*             0 - failure                                                    */
/* History    :                                                               */
/*             31/08/98 SMJS  Initial Implementation                          */
/******************************************************************************/
int EcoString_changeStr(ECOSTRTABLE *EcoSTabP, ECOSTRING *To, char *From,int StartPos) {
  char *TmpStr;

  if (!EcoString_subtractOne(EcoSTabP,*To)) {
    Error_trace("EcoString_changeStr",NULL);
    return 0;
  }
  *To = NULL;
  /* Generate the new string */
  if (!StrUtil_copyString(&TmpStr,From,StartPos)) {
    Error_trace("EcoString_changeStr",NULL);
    return 0;
  }
  if (!EcoString_getPointer(EcoSTabP,To,TmpStr)) {
    Error_trace("EcoString_changeStr",NULL);
    return 0;
  }

/* Return success */
  return 1;
}

/******************************************************************************/
/* Routine    :                                                               */
/*             EcoString_chkCheckSum()                                        */
/* Role       :                                                               */
/*             Check a very simple checksum for a string.                     */
/* Arguments  :                                                               */
/*             String - The string to generate checksum for                   */
/* Returns    :                                                               */
/*             1 - Checksum OK.                                               */
/*             0 - Checksum NOT OK.                                           */
/* History    :                                                               */
/*             31/08/98 SMJS  Initial Implementation                          */
/******************************************************************************/
int EcoString_chkCheckSum(ECOSTRING String,int CheckSum) {
  int ChkSum = 0;
  int i;
  int Mult;
  int Len;

  if (String == NULL)
    return CheckSum==0;

  Len = strlen(String);
  for (i=0,Mult=1;i<Len;i++, Mult=(-Mult)) {
    ChkSum += (String[i]*Mult);
  }
  return ChkSum==CheckSum;
}

/******************************************************************************/
/* Routine    :                                                               */
/*             EcoString_copyNStr()                                           */
/* Role       :                                                               */
/*             Get a pointer to a (possibly new) EcoString                    */
/* Arguments  :                                                               */
/*             EcoSTabP - Pointer to ECOSTRTABLE to add to                    */
/*             To       - Pointer to ECOSTRING (RETURNED)                     */
/*             From     - String to copy from                                 */
/*             StartPos - the first position in From to copy                  */
/*             Length   - the maximum length to get                           */
/* Returns    :                                                               */
/*             1 - success                                                    */
/*             0 - failure                                                    */
/* History    :                                                               */
/*             31/08/98 SMJS  Initial Implementation                          */
/******************************************************************************/
int EcoString_copyNStr(ECOSTRTABLE *EcoSTabP, ECOSTRING *To, char *From,int StartPos,
                       int Length) {
  char *TmpStr;

  *To = NULL;
  /* Generate the string */
  if (!StrUtil_copyNString(&TmpStr,From,StartPos,Length)) {
    Error_trace("EcoString_copyNStr",NULL);
    return 0;
  }
  if (!EcoString_getPointer(EcoSTabP,To,TmpStr)) {
    Error_trace("EcoString_copyNStr",NULL);
    return 0;
  }

/* Return success */
  return 1;
}

/******************************************************************************/
/* Routine    :                                                               */
/*             EcoString_copyNTStr()                                          */
/* Role       :                                                               */
/*             Get a pointer to a (possibly new) EcoString                    */
/* Arguments  :                                                               */
/*             EcoSTabP - Pointer to ECOSTRTABLE to add to                    */
/*             To       - Pointer to ECOSTRING (RETURNED)                     */
/*             From     - String to copy from                                 */
/*             StartPos - the first position in From to copy                  */
/*             Length   - the maximum length to get                           */
/* Returns    :                                                               */
/*             1 - success                                                    */
/*             0 - failure                                                    */
/* History    :                                                               */
/*             31/08/98 SMJS  Initial Implementation                          */
/******************************************************************************/
int EcoString_copyNTStr(ECOSTRTABLE *EcoSTabP, ECOSTRING *To, char *From,int StartPos,
                        int Length) {
  char *TmpStr;

  *To = NULL;

/* Generate the string */
  if (!StrUtil_copyNTString(&TmpStr,From,StartPos,Length)) {
    Error_trace("EcoString_copyNTStr",NULL);
    return 0;
  }
  if (!EcoString_getPointer(EcoSTabP,To,TmpStr)) {
    Error_trace("EcoString_copyNTStr",NULL);
    return 0;
  }

/* Return success */
  return 1;
}

/******************************************************************************/
/* Routine    :                                                               */
/*             EcoString_copyStr()                                            */
/* Role       :                                                               */
/*             Get a pointer to a (possibly new) EcoString                    */
/* Arguments  :                                                               */
/*             EcoSTabP - Pointer to ECOSTRTABLE to add to                    */
/*             To       - Pointer to ECOSTRING (RETURNED)                     */
/*             From     - String to copy from                                 */
/*             StartPos - the first position in From to copy                  */
/* Returns    :                                                               */
/*             1 - success                                                    */
/*             0 - failure                                                    */
/* History    :                                                               */
/*             31/08/98 SMJS  Initial Implementation                          */
/******************************************************************************/
int EcoString_copyStr(ECOSTRTABLE *EcoSTabP, ECOSTRING *To, char *From,int StartPos) {
  char *TmpStr;

  *To = NULL;

/* Generate the string */
  if (!StrUtil_copyString(&TmpStr,From,StartPos)) {
    Error_trace("EcoString_copyStr",NULL);
    return 0;
  }
  if (!EcoString_getPointer(EcoSTabP,To,TmpStr)) {
    Error_trace("EcoString_copyStr",NULL);
    return 0;
  }

/* Return success */
  return 1;
}

/******************************************************************************/
/* Routine    :                                                               */
/*             EcoString_copyTStr()                                           */
/* Role       :                                                               */
/*             Get a pointer to a (possibly new) EcoString                    */
/* Arguments  :                                                               */
/*             EcoSTabP - Pointer to ECOSTRTABLE to add to                    */
/*             To       - Pointer to ECOSTRING (RETURNED)                     */
/*             From     - String to copy from                                 */
/*             StartPos - the first position in From to copy                  */
/* Returns    :                                                               */
/*             1 - success                                                    */
/*             0 - failure                                                    */
/* History    :                                                               */
/*             31/08/98 SMJS  Initial Implementation                          */
/******************************************************************************/
int EcoString_copyTStr(ECOSTRTABLE *EcoSTabP, ECOSTRING *To, char *From,int StartPos) {
  char *TmpStr;

  *To = NULL;
/* Generate the string */
  if (!StrUtil_copyTString(&TmpStr,From,StartPos)) {
    Error_trace("EcoString_copyTStr",NULL);
    return 0;
  }
  if (!EcoString_getPointer(EcoSTabP,To,TmpStr)) {
    Error_trace("EcoString_copyTStr",NULL);
    return 0;
  }

/* Return success */
  return 1;
}

/******************************************************************************/
/* Routine    :                                                               */
/*             EcoString_delStr()                                             */
/* Role       :                                                               */
/*             Remove an ECOSTRING from an ECOSTRTABLE                        */
/* Arguments  :                                                               */
/*             EcoSTabP - Pointer to ECOSTRTABLE                              */
/*             String   - ECOSTRING                                           */
/*             StrInd   - Index of String in EcoSTabP->CHashTab->Strings      */
/* Returns    :                                                               */
/*             1 - success                                                    */
/*             0 - failure                                                    */
/* History    :                                                               */
/*             31/08/98 SMJS  Initial Implementation                          */
/******************************************************************************/
int EcoString_delStr(ECOSTRTABLE *EcoSTabP, ECOSTRING String,int StrInd) {
  int LetInd;

  if (EcoSTabP->UseCount[StrInd]>0) {
    Error_write(EECOSTR,"EcoString_delStr",ERR_SEVERE,
             "Trying to delete ECOSTRING %s which is in use %d times",
             String, EcoSTabP->UseCount[StrInd]);
    return 0;
  }

/* Look for it in array */
  if (!EcoSTabP->CHashTab->HashFunc(String,&LetInd)) {
    Error_trace("EcoString_delStr",NULL);
    return 0;
  } else {
    EcoSTabP->CheckSums[StrInd]=0;
    if (!CHash_removeStr(EcoSTabP->CHashTab,String,StrInd,LetInd)) {
      Error_trace("EcoString_delStr",NULL);
      return 0;
    }
  }

/* Return success */
  return 1;
}

/******************************************************************************/
/* Routine    :                                                               */
/*             EcoString_freeStr()                                            */
/* Role       :                                                               */
/*             Reduce the UseCount and possibly free an ECOSTRING. This       */
/*             routine is just a wrapper for EcoString_subtractOne at moment  */
/* Arguments  :                                                               */
/*             EcoSTabP - Pointer to ECOSTRTABLE                              */
/*             String   - ECOSTRING                                           */
/* Returns    :                                                               */
/*             1 - success                                                    */
/*             0 - failure                                                    */
/* History    :                                                               */
/*             31/08/98 SMJS  Initial Implementation                          */
/******************************************************************************/
int EcoString_freeStr(ECOSTRTABLE *EcoSTabP, ECOSTRING String) {
  if (!EcoString_subtractOne(EcoSTabP,String)) {
    Error_trace("EcoString_freeStr",NULL);
    return 0;
  }

/* Return success */
  return 1;
}

/******************************************************************************/
/* Routine    :                                                               */
/*             EcoString_genCheckSum()                                        */
/* Role       :                                                               */
/*             Generate a very simple checksum for a string.                  */
/* Arguments  :                                                               */
/*             String - The string to generate checksum for                   */
/* Returns    :                                                               */
/*             the checksum (ALWAYS)                                          */
/* History    :                                                               */
/*             31/08/98 SMJS  Initial Implementation                          */
/******************************************************************************/
int EcoString_genCheckSum(char *String) {
  int ChkSum = 0;
  int i;
  int Mult;
  int Len;

  if (String == NULL) return 0;

  Len = strlen(String);
  for (i=0,Mult=1;i<Len;i++, Mult=(-Mult)) {
    ChkSum += (String[i]*Mult);
  }
/* Return checksum */
  return ChkSum;
}

/******************************************************************************/
/* Routine    :                                                               */
/*             EcoString_getInfo()                                            */
/* Role       :                                                               */
/*             Display information about an ECOSTRTABLE                       */
/* Arguments  :                                                               */
/*             EcoSTabP - Pointer to ECOSTRTABLE to display                   */
/* Returns    :                                                               */
/*             1 - success                                                    */
/*             0 - failure                                                    */
/* History    :                                                               */
/*             31/08/98 SMJS  Initial Implementation                          */
/******************************************************************************/
int EcoString_getInfo(ECOSTRTABLE *EcoSTabP) {
  int i;
  int TotalUsed = 0;
  int TotalSaved = 0;
  int Len;
  int nElem = 0;

  for (i=0;i<EcoSTabP->CHashTab->NElement;i++) {
    if (EcoSTabP->CHashTab->Strings[i]) {
      nElem++;
      Len=strlen(EcoSTabP->CHashTab->Strings[i])+1;
      TotalUsed+=Len;
      if (EcoSTabP->UseCount[i]>1) {
        TotalSaved+=Len*(EcoSTabP->UseCount[i]-1) -
                    (sizeof(CHASHARRAY)+sizeof(int)+sizeof(int));
      } else {
        TotalSaved-=(sizeof(CHASHARRAY)+sizeof(int)+sizeof(int));
      }
    } else {
      TotalSaved-=(sizeof(char *)+sizeof(int)+sizeof(int));
    }
  }
  // Note NElement is currently not necessarily number of strings because of removals printf("Number of strings in EcoS Table:  %10d\n",EcoSTabP->CHashTab->NElement);
  printf("Number of strings in EcoS Table:  %10d\n", nElem);
  printf("Total bytes used to hold strings: %10d\n", TotalUsed);
  printf("Total bytes SAVED in using table: %10d\n", TotalSaved);

/* Return success */
  return 1;
}

/******************************************************************************/
/* Routine    :                                                               */
/*             EcoString_getPointer()                                         */
/* Role       :                                                               */
/*             Get a pointer to a (possibly new) EcoString                    */
/* Arguments  :                                                               */
/*             EcoSTabP - Pointer to ECOSTRTABLE to add to                    */
/*             To       - Pointer to ECOSTRING (RETURNED)                     */
/*             From     - String to copy from                                 */
/* Returns    :                                                               */
/*             1 - success                                                    */
/*             0 - failure                                                    */
/* History    :                                                               */
/*             31/08/98 SMJS  Initial Implementation                          */
/******************************************************************************/
int EcoString_getPointer(ECOSTRTABLE *EcoSTabP, ECOSTRING *To, char *From) {
  int StrInd;

/* Look for it in array */
  if (!CHash_find(From,EcoSTabP->CHashTab,&StrInd)) {
/* If not found add it to array */
    if (!EcoString_addStr(EcoSTabP,From,&StrInd)) {
      Error_trace("EcoString_getPointer",NULL);
      free(From);
      return 0;
    }
  } else {
    //printf("Old string %s\n",From);
    (EcoSTabP->UseCount[StrInd])++;
    free(From);
  }
/* Set the pointer */
  *To = EcoSTabP->CHashTab->Strings[StrInd];

/* Return success */
  return 1;
}

/******************************************************************************/
/* Routine    :                                                               */
/*             EcoString_initTable()                                          */
/* Role       :                                                               */
/*             Initialise an ECOSTRTABLE                                      */
/* Arguments  :                                                               */
/*             EcoSTabP - Pointer to ECOSTRTABLE to initialise (RETURNED)     */
/* Returns    :                                                               */
/*             1 - success                                                    */
/*             0 - failure                                                    */
/* History    :                                                               */
/*             31/08/98 SMJS  Initial Implementation                          */
/******************************************************************************/
int EcoString_initTable(ECOSTRTABLE **EcoSTabP) {
  if (((*EcoSTabP) = (ECOSTRTABLE *)calloc(1,sizeof(ECOSTRTABLE)))==NULL) {
    Error_write(ECALLERR,"EcoString_init",ERR_SEVERE,"*EcoSTabP");
    return 0;
  }
  if (!CHash_setup(NULL,0,&((*EcoSTabP)->CHashTab),CHash_fourLets,1024)) {
    Error_trace("EcoString_init",NULL);
    return 0;
  }
  return 1;
}

/******************************************************************************/
/* Routine    :                                                               */
/*             EcoString_subtractOne()                                        */
/* Role       :                                                               */
/*             Subtract one from the UseCount for an ECOSTRING                */
/* Arguments  :                                                               */
/*             EcoSTabP - Pointer to ECOSTRTABLE                              */
/*             String   - ECOSTRING                                           */
/* Returns    :                                                               */
/*             1 - success                                                    */
/*             0 - failure                                                    */
/* History    :                                                               */
/*             31/08/98 SMJS  Initial Implementation                          */
/******************************************************************************/
int EcoString_subtractOne(ECOSTRTABLE *EcoSTabP, ECOSTRING String) {
  int StrInd;

/* Look for it in array */
  if (!CHash_find(String,EcoSTabP->CHashTab,&StrInd)) {
    Error_write(EECOSTR,"EcoString_subtractOne",ERR_SEVERE,
             "%s is not in EcoString table",String);
    return 0;
  } else {
    if (EcoSTabP->UseCount[StrInd]<=0) {
      Error_write(EECOSTR,"EcoString_subtractOne",ERR_SEVERE,
               "Use count for %s is invalid (= %d). StrInd = %d",
               String,EcoSTabP->UseCount[StrInd],StrInd);
      return 0;
    }
    (EcoSTabP->UseCount[StrInd])--;
    if (EcoSTabP->UseCount[StrInd] == 0) {
      if (!EcoString_delStr(EcoSTabP,String,StrInd)) {
        Error_trace("EcoString_subtractOne",NULL);
        return 0;
      }
    }
  }

/* Return success */
  return 1;
}
