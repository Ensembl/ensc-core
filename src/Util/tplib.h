/*
    TPatterns
  
    Copyright (C) T J R Cutts 1998-2002

    timc@chiark.greenend.org.uk

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <limits.h>

/* Translation matrix */
typedef char txmatrix[4][4][4];
extern txmatrix matrix;

/* Functions in translate.c */
void make_revmatrix(void);
void initbasebits(void);
int compilemx(char *);
void translate(char *, char **, int*);
#ifdef NVSN
void rev_comp(char *, char *);
#endif

/* Error stuff in error.c */
typedef enum {
  TPERR_SUCCESS = 0,
  TPERR_MEM,
  TPERR_MXEOF,
  TPERR_REDIR
} tperror;

void tp_error(tperror);
