#ifndef __SEQEDIT_H__
#define __SEQEDIT_H__

#include <stdio.h>
#include "DataModelTypes.h"
#include "EnsRoot.h"
#include "EcoString.h"

OBJECTFUNC_TYPES(SeqEdit)

typedef struct SeqEditFuncsStruct {
  OBJECTFUNCS_DATA(SeqEdit)
} SeqEditFuncs;

#define FUNCSTRUCTTYPE SeqEditFuncs
struct SeqEditStruct {
  ENSROOT_DATA
  char *altSeq;
  ECOSTRING name;
  ECOSTRING code;
  ECOSTRING description;
  long start;
  long end;
};
#undef FUNCSTRUCTTYPE

SeqEdit *SeqEdit_new(long start, long end, char *altSeq, char *name, char *desc, char *code);
SeqEdit *SeqEdit_newFromAttribute(Attribute *attrib);

long SeqEdit_getLengthDiff(SeqEdit *seqEd);
Attribute *SeqEdit_getAttribute(SeqEdit *seqEd);
char *SeqEdit_applyEdit(SeqEdit *seqEd, char *seq);

int SeqEdit_reverseStartCompFunc(const void *one, const void *two);

ECOSTRING SeqEdit_setName(SeqEdit *seqEdit, char *name);
#define SeqEdit_getName(seqEd) (seqEd)->name

ECOSTRING SeqEdit_setCode(SeqEdit *seqEdit, char *code);
#define SeqEdit_getCode(seqEd) (seqEd)->code

ECOSTRING SeqEdit_setDescription(SeqEdit *seqEdit, char *description);
#define SeqEdit_getDescription(seqEd) (seqEd)->description

ECOSTRING SeqEdit_setAltSeq(SeqEdit *seqEdit, char *altSeq);
#define SeqEdit_getAltSeq(seqEd) (seqEd)->altSeq

long SeqEdit_setStart(SeqEdit *seqEdit, long start);
#define SeqEdit_getStart(seqEd) (seqEd)->start

long SeqEdit_setEnd(SeqEdit *seqEdit, long end);
#define SeqEdit_getEnd(seqEd) (seqEd)->end

void SeqEdit_free(SeqEdit *seqEdit);

#ifdef __SEQEDIT_MAIN__
  SeqEditFuncs
    seqEditFuncs = {
                    SeqEdit_free,
                    NULL, // shallowCopy
                    NULL  // deepCopy
                   };
#else
  extern SeqEditFuncs seqEditFuncs;
#endif



#endif
