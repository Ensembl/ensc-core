#include "SliceAdaptor.h"
#include "BaseAdaptor.h"
#include "MysqlUtil.h"


SliceAdaptor *SliceAdaptor_new(DBAdaptor *dba) {
  SliceAdaptor *sa;

  if ((sa = (SliceAdaptor *)calloc(1,sizeof(SliceAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for SliceAdaptor\n");
    return NULL;
  }
  BaseAdaptor_init((BaseAdaptor *)sa, dba, SLICE_ADAPTOR);

  return sa;
}

Slice *SliceAdaptor_fetchByChrStartEnd(SliceAdaptor *sa, char *chr, int start, int end) {
  Slice *slice;
  char *assemblyType;

  if (!chr) {
    fprintf(stderr,"ERROR: chromosome name argument must be defined and not ''\n");
    exit(1);
  }

  if(start > end) {
    fprintf(stderr,"ERROR: start must be less than end: parameters %s:%d:%d\n",chr,start,end);
    exit(1);
  }

  assemblyType = DBAdaptor_getAssemblyType(sa->dba);

  slice = Slice_new(chr, start, end, 1, assemblyType, sa, 0, FALSE);

  return slice;
}

