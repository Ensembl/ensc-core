#ifndef __ECOSTR_H__
#define __ECOSTR_H__

#include "CHash.h"

typedef char *                 ECOSTRING;
typedef struct EcoStrTabStruct ECOSTRTABLE;

struct EcoStrTabStruct {
  CHASHTABLE *CHashTab;
  int        *CheckSums;
  int        *UseCount;
};

#define EcoString_strcmp(X,Y) !((X) == (Y))

/* prototypes */
int EcoString_addStr(ECOSTRTABLE *EcoSTabP, char *String, int *StrInd);
int EcoString_changeStr(ECOSTRTABLE *EcoSTabP,ECOSTRING *To,char *From,int StartPos);
int EcoString_chkCheckSum(ECOSTRING String,int CheckSum);
int EcoString_copyStr(ECOSTRTABLE *EcoSTabP,ECOSTRING *To,char *From,int StartPos);
int EcoString_copyTStr(ECOSTRTABLE *EcoSTabP,ECOSTRING *To,char *From,int StartPos);
int EcoString_copyNStr(ECOSTRTABLE *EcoSTabP,ECOSTRING *To,char *From,int StartPos,
		       int Length);
int EcoString_copyNTStr(ECOSTRTABLE *EcoSTabP,ECOSTRING *To,char *From,int StartPos,
		        int Length);
int EcoString_delStr(ECOSTRTABLE *EcoSTabP,ECOSTRING String,int StrInd);
int EcoString_freeStr(ECOSTRTABLE *EcoSTabP, ECOSTRING String);
int EcoString_genCheckSum(char *String);
int EcoString_getInfo(ECOSTRTABLE *EcoSTabP);
int EcoString_getPointer(ECOSTRTABLE *EcoSTabP, ECOSTRING *To, char *From);
int EcoString_initTable(ECOSTRTABLE **EcoSTabP);
int EcoString_subtractOne(ECOSTRTABLE *EcoSTabP, ECOSTRING String);

#endif /* __ECOSTR_H__ */
