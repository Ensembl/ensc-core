#include "Clone.h"
#include "CloneAdaptor.h"

Clone *Clone_new() {
  Clone *cl;

  if ((cl = (Clone *)calloc(1,sizeof(Clone))) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for cl\n");
    return NULL;
  }

  return cl;
}

char *Clone_setName(Clone *cl, char *name) {
  if ((cl->name = (char *)malloc(strlen(name)+1)) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for clone name\n");
    return NULL;
  }

  strcpy(cl->name,name);

  return cl->name;
}

char *Clone_setEmblAcc(Clone *cl, char *acc) {
  if ((cl->emblAcc = (char *)malloc(strlen(acc)+1)) == NULL) {
    fprintf(stderr,"ERROR: Failed allocating space for clone acc\n");
    return NULL;
  }

  strcpy(cl->emblAcc,acc);

  return cl->emblAcc;
}

