#ifndef __CIGARSTRUTIL_H__
#define __CIGARSTRUTIL_H__

#include "Vector.h"

char *CigarStrUtil_reverse(char *oldCigarString, int len);
Vector *CigarStrUtil_getPieces(char *cigarString);

#endif
