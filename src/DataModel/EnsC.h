#ifndef __ENSC_H__
#define __ENSC_H__

#include <stdio.h>

#ifdef TRUE
#undef TRUE
#endif
#define TRUE 1

#ifdef FALSE
#undef FALSE
#endif
#define FALSE 0

typedef enum CopyDepthEnum {
  NONE_DEPTH,
  SHALLOW_DEPTH,
  DEEP_DEPTH
} CopyDepth;

#endif
