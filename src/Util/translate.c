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

#include "tplib.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/* The following matrix encodes the Universal genetic code.  During
   execution any matrix supplied via the -m option will overwrite it */

txmatrix matrix;

txmatrix standard_matrix = {
  { "KNKN", "TTTT", "RSRS", "IIMI" },
  { "QHQH", "PPPP", "RRRR", "LLLL" },
  { "EDED", "AAAA", "GGGG", "VVVV" },
  { "*Y*Y", "SSSS", "*CWC", "LFLF" }
};

txmatrix ct2id_matrix = {
  { "KNKN", "TTTT", "*S*S", "MIMI" },
  { "QHQH", "PPPP", "RRRR", "LLLL" },
  { "EDED", "AAAA", "GGGG", "VVVV" },
  { "*Y*Y", "SSSS", "WCWC", "LFLF" }
};

txmatrix revmatrix;

static int initTable = 0;

/* make_revmatrix() creates the reverse complemented matrix by
   swapping bits 4 and 5 with bits 0 and 1, inverting the bits, and
   then keeping only bits 0-5 */

void make_revmatrix(void) {
  int n, r;

  for (n = 0; n<64; n++) {
    r = (~((n & 0x0C) |
          ((n & 0x03) << 4) |
          ((n >> 4) & 0x03))) & 0x3F;

    ((char *)revmatrix)[r] = ((char*)matrix)[n];
  }

}

/* compilemx() reads a file of codon -> amino acid translations and
   places them in matrix, above, for later use by the translate()
   function */

int compilemx(char *filename) {
  char buf[80];
  char path[PATH_MAX];

  FILE *f;
  int a, b, c, res;
  char ta,tb,tc,tx,*r;

  if ((f=fopen(filename, "r"))==NULL) {
    /* Don't bother with the environment if an absolute path was
       given */
    if (filename[0] != '/') {
      if ((r = getenv("ENSC_DATADIR"))!=NULL) {
        sprintf(path, "%s/%s", r, filename);
        f=fopen(path, "r");
      }
/*
      if (f == NULL) {
        sprintf(path, "%s/%s", TFLIB, filename);
        f=fopen(path, "r");
      }
*/
    }
  }

  if (f == NULL) {
    perror("Could not open translation matrix file:");
    exit(1);
  }

  for (a = 0; a<4; a++) {
    for (b = 0; b<4; b++) {
      for (c = 0; c<4; c++) {
        do {
          r = fgets(buf, sizeof(buf), f);
          if (r==NULL) {
            /* Premature end of file encountered */
            fprintf(stderr, "Premature end of matrix file encountered for %s\n",filename);
            exit(1);
          }
          res = sscanf(buf, "%c%c%c %c", &ta, &tb, &tc, &tx);
        } while (res!=4);

        matrix[a][b][c]=tx;
      }
    }
  }

  fclose(f);
  return 0;
}

/* initbasebits() initialises the lookup table for the translate()
   function */

int basebits[256];

int comphash[256];

void initbasebits(void)
{

  /* First of all the hash for N->P translation */

  memset(basebits, 255, sizeof(int)*256);
  basebits['A'] = 0;
  basebits['C'] = 1;
  basebits['G'] = 2;
  basebits['T'] = 3;
  basebits['U'] = 3;
  basebits['a'] = 0;
  basebits['c'] = 1;
  basebits['g'] = 2;
  basebits['t'] = 3;
  basebits['u'] = 3;

  /* And now the hash for reverse complementation */
  memset(comphash, 'x', sizeof(int)*256);
  comphash['A'] = 'T';
  comphash['B'] = 'V';
  comphash['C'] = 'G';
  comphash['D'] = 'H';
  comphash['G'] = 'C';
  comphash['H'] = 'D';
  comphash['K'] = 'M';
  comphash['M'] = 'K';
  comphash['N'] = 'N';
  comphash['R'] = 'Y';
  comphash['S'] = 'S';
  comphash['T'] = 'A';
  comphash['U'] = 'A';
  comphash['V'] = 'B';
  comphash['X'] = 'X';
  comphash['Y'] = 'Y';
  comphash['.'] = '.';
  comphash['~'] = '~';

}

/* translate() a nucleic acid sequence in all six reading frames
   simultaneously.  out should be an array of six (char *) pointers,
   all of which must be large enough to hold the translated sequence.
   No sanity checking is performed here. l is an array of six
   integers, which will take the lengths of the translated sequences.

   Thanks to Peter Benie for help optimising the code */

void translate(char *in, char **out, int *l, int codonTableId) {
  int n, rem, index;
  char *p, *end, *r3, *r4, *r5;
  char *r0 = out[0];
  char *r1 = out[1];
  char *r2 = out[2];

  if (initTable != codonTableId) {
 // Hack for now
    switch (codonTableId) {
      case 1:
        memcpy(matrix, standard_matrix, sizeof(txmatrix));
        break;
      case 2:
        memcpy(matrix, ct2id_matrix, sizeof(txmatrix));
        break;
      default:
        fprintf(stderr,"Error: Currently unsupported codon table id (%d)\n", codonTableId);
        exit(1);
    }

    initbasebits();
    make_revmatrix();
    initTable = codonTableId;
  }

  n = strlen(in);
  end = &in[n]-5;

  /* How much of a part codon do we have at the end */
  rem = n % 3;

  /* How many residues long will the frame 0 sequence be? */
  n /= 3;

  l[0] = l[1] = l[2] = l[3] = l[4] = l[5] = n;

  /* Position the reverse reading frame pointers where the
     protein sequences will end */

  r3 = &out[3][n];
  r4 = &out[4][n];
  r5 = &out[5][n];

  /* If there aren't two bases at the end, the translations
     will not all be the same length; one or two of them
     will be one residue shorter */

  if (rem != 2) {
    r5--;
    l[5]--;
    l[2]--;
    if (rem == 0) {
      r4--;
      l[4]--;
      l[1]--;
    }
  }

  *r3 = *r4 = *r5 = '\0';
  r3--; r4--; r5--;

  index = (basebits[(int)in[0]] << 2) | basebits[(int)in[1]];

  for (p = in; p<=end; p+=3) {
    /* Frame 0/3 */

    index = ((index << 2) & 0x3C) | basebits[(int)p[2]];

// Hack for N's
    if (p[0] == 'N' || p[1] == 'N' || p[2] == 'N') index = index | 0x80;

    if (index & 0x80) {
      *r0++ = 'X';
      *r3-- = 'X';
    } else {
      *r0++ = ((char *)matrix)[index];
      *r3-- = ((char *)revmatrix)[index];
    }
    /* Frame 1/4 */

    index = ((index << 2) & 0x3C) | basebits[(int)p[3]];

    if (index & 0x80) {
      *r1++ = 'X';
      *r4-- = 'X';
    } else {
      *r1++ = ((char *)matrix)[index];
      *r4-- = ((char *)revmatrix)[index];
    }

    /* Frame 2/5 */

    index = ((index << 2) & 0x3C) | basebits[(int)p[4]];

    if (index & 0x80) {
      *r2++ = 'X';
      *r5-- = 'X';
    } else {
      *r2++ = ((char *)matrix)[index];
      *r5-- = ((char *)revmatrix)[index];
    }
  }

  /* We may, at this point, have one or two codons still untranslated,
     so we need to deal with those */

  if (rem<2) {
    index = ((index << 2) & 0x3C) | basebits[(int)p[2]];

    if (index & 0x80) {
      *r0++ = 'X';
      *r3-- = 'X';
    } else {
      *r0++ = ((char *)matrix)[index];
      *r3-- = ((char *)revmatrix)[index];
    }

    if (rem==1) {
      index = ((index << 2) & 0x3C) | basebits[(int)p[3]];

      if (index & 0x80) {
        *r1++ = 'X';
        *r4-- = 'X';
      } else {
        *r1++ = ((char *)matrix)[index];
        *r4-- = ((char *)revmatrix)[index];
      }
    }
  }
    
  *r0 = *r1 = *r2 = '\0';
}

/* rev_comp() produces the reverse complement of a nucleic acid
   sequence.  It is only used when comparing a nucleotide pattern
   against a nucleotide db */

void rev_comp(char *in, char *out, int length) {
  register char *p = in;
  register char *q;
  
  q = &out[length];
  *q-- ='\0';
  
  while (q>=out)
    *q-- = comphash[*p++];
}
