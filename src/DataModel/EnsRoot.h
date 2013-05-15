#ifndef __ENSROOT_H__
#define __ENSROOT_H__

#include "Object.h"

#define ENSROOT_DATA OBJECT_DATA

#define EnsRoot_free(ef) Object_free((ef))

#define EnsRoot_shallowCopy(sf) Object_shallowCopy((sf))
#define EnsRoot_deepCopy(sf) Object_deepCopy((sf))

#endif
