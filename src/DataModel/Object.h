#ifndef __OBJECT_H__
#define __OBJECT_H__

#include "Class.h"

typedef struct ObjectStruct Object;

#define OBJECT_DATA \
  ClassType objectType;

struct ObjectStruct {
  OBJECT_DATA
};

#endif
