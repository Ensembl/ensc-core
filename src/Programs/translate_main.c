/*
 * Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
 * Copyright [2016-2017] EMBL-European Bioinformatics Institute
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "tplib.h"

#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define LINELENGTH 60
#define STOPCODON "*"
#define MIN_LENGTH 20
#define NUMSEQ 20
#define TRUE 1
#define FALSE 0

#define FLAG_KEEPSTOPS 2
#define FLAG_MET 4
#define FLAG_BESTONLY 8

#define OPTIONS "ahl:mb"

int minimum_len = MIN_LENGTH;

static char usage[] = "\
Usage: translate [-options] <seqfile>\n\
   Translate a nucleic acid sequence into protein ORFs.\n\
   Available options are:\n\
   -a            : translate in full, with stops; no individual ORFs\n\
   -l <minlen>   : report only ORFs greater than minlen (default 20)\n\
   -m            : require ORFs to start with AUG/Met\n\
   -b            : only output longest open reading frame\n\
   -h            : help\n";

typedef struct {
  char name[1024];
  char *seq;
  int length;
} SEQ;

void print_sequence (char *header, char *s, int n) {
  int i = LINELENGTH;
  printf("%s\n", header);
  while (i < n) {
    if (fwrite(s, sizeof(char), LINELENGTH, stdout) == 0) exit(EXIT_FAILURE);
    if (fputc('\n', stdout) == EOF) exit(EXIT_FAILURE);
    i += LINELENGTH;
    s += LINELENGTH;
  }
  printf("%s\n", s);
}

int seqComp(const void *a, const void *b) {
  SEQ **seq1 = (SEQ **)a;
  SEQ **seq2 = (SEQ **)b;

  if ((*seq1)->length > (*seq2)->length) {
    return -1;
  } else if ((*seq1)->length < (*seq2)->length) {
    return 1;
  } else {
    return 0;
  }
}

void process_sequence (char *s, char **frm, int *lengths, char *title, int flags) {
  int n;
  char *orf;
  SEQ *seq;
  SEQ **seqs = NULL;
  int nseq = 0;
  int seq_length = strlen(s);
  int orfnumber = 1;
  int start;
  int end;
  int orf_length;
  char header[1024];

  translate(s, frm, lengths, 1, seq_length);
  for (n = 0; n<6; n++){
    if (flags & FLAG_KEEPSTOPS) {
      if (n < 3) {
        start = n+1;
        end = lengths[n]*3+start-1;
      }
      else {
        start = seq_length-n+3;
        end = start-lengths[n]*3+1;
      }
      sprintf(header, "%s.%d\tlength %i, nt %i..%i", title, n, lengths[n], start, end);
      print_sequence(header, frm[n], lengths[n]);
    }
    else {
      orf = strtok(frm[n], STOPCODON);
      while (orf != NULL && *orf != '\0') {
        if (flags & FLAG_MET) {
          while (*orf != 'M' && *orf != '\0') orf++;
        }
        if (*orf != '\0') {
          orf_length = strlen(orf);
          if (orf_length > minimum_len) {
            if (nseq == 0) {
              if ((seqs = calloc(NUMSEQ,sizeof(SEQ *))) == NULL) {
                exit(EXIT_FAILURE);
              }
            }
            else if (nseq%NUMSEQ == 0) {
              seqs = realloc(seqs,(nseq+NUMSEQ)*sizeof(SEQ *));
            }
            seqs[nseq++] = calloc(1,sizeof(SEQ));
            seq = seqs[nseq-1];

            start = (orf - frm[n]) * 3 + 1;
            if (n < 3) {
              start += n;
              end = start+orf_length*3-1;
            } else {
              start = seq_length-start-n+4;
              end = start-orf_length*3+1;
            }

            seq->length = orf_length;
            sprintf(seq->name,"%s.%d\tlength %i, nt %i..%i", title, orfnumber, orf_length, start, end);
            seq->seq = calloc(orf_length+1,sizeof(char));
            strcpy(seq->seq,orf);

            orfnumber++;
          }
        }
        orf = strtok(NULL, STOPCODON);
      }
    }
  }
  *s = '\0';
  if (nseq) {
    qsort(seqs, nseq, sizeof(SEQ*), seqComp);
    if (flags & FLAG_BESTONLY) {
      print_sequence(seqs[0]->name, seqs[0]->seq, seqs[0]->length);
    }
    else {
      for (n = 0; n < nseq; n++) {
        print_sequence(seqs[n]->name, seqs[n]->seq, seqs[n]->length);
      }
    }
    for (n = 0; n < nseq; n++) {
      free(seqs[n]->seq);
      free(seqs[n]);
    }
    free(seqs);
  }
}

void process_file(char *filename, int flags) {
  char *s;
  char *p;
  char *frm[6];
  static char buf[1024];
  static char title[1024];
  int sz = 8192;
  int seq_length;
  int lengths[6];
  int i;
  int n;
  FILE *f;

  s = (char *)malloc(sz);
  if (!s) {
    exit(EXIT_FAILURE);
  }
  for ( i = 0; i < 6; i++) {
    frm[i] = (char *)malloc(sz);

    if (!frm[i]) {
      exit(EXIT_FAILURE);
    }
  }

  *s = '\0';
  seq_length = 0;

  f = fopen(filename, "r");
  if (f == NULL) exit(EXIT_FAILURE);

  while (!feof(f)) {
    p = fgets(buf, sizeof(buf), f);
    if (p) {
      if (*p == '>') {
        /* We have found a comment line */
        if (*s != '\0') {
          s[seq_length] = '\0';
          /* We now have a complete sequence, so
             process it */
          process_sequence(s, frm, lengths, title, flags);
          seq_length = 0;
        }
        buf[strlen(buf)-1] = '\0';
        strcpy(title, buf);
      }
      else {
        /* This is a sequence line, so add it to s */
        n = strlen(p)-1;
        if ((seq_length+n+3)>sz) {
          sz <<= 1;
          s = (char *)realloc(s, sz);
          if (!s) {
            exit(EXIT_FAILURE);
          }
          for (i = 0; i < 6; i++) {
            frm[i] = (char *)realloc(frm[i], sz);
            if (!frm[i]) {
              exit(EXIT_FAILURE);
            }
          }
        }
        memcpy(&s[seq_length], buf, n);
        seq_length += n;
      }
    }
  }

  if (*s != '\0') {
    s[seq_length] = '\0';
    /* We now have a complete sequence, so
       process it */
    process_sequence(s, frm, lengths, title, flags);
    seq_length = 0;
  }

  if (filename)
    fclose(f);

  for (i = 0; i < 6; i++) {
    free(frm[i]);
  }
  free(s);
}

int main(int argc, char *argv[]) {
  int n;
  int optchar;
  int flags = 0; //Flag is set to 0 to get any ORF which is bigger than MIN_LENGTH
  int nopts = 1;
  char nameStr[1024];

  while ((optchar = getopt(argc, argv, OPTIONS)) != -1) {
    switch (optchar) {
      case 'a': flags += FLAG_KEEPSTOPS; break;
      case 'm': flags += FLAG_MET; break;
      case 'b': flags += FLAG_BESTONLY; break;
      case 'l': minimum_len = atoi(optarg); ++nopts; break;

      case 'h':
                fprintf(stderr, "%s\n", usage);
                exit(EXIT_SUCCESS);
      default:
                fprintf(stderr, "%s\n", usage);
                exit(EXIT_FAILURE);
    }
    ++nopts;
  }
  if (nopts == argc) {
    fprintf(stderr, "%s\n", usage);
    exit(EXIT_FAILURE);
  }
  make_revmatrix();
  for (n = nopts; n < argc; n++) {
    process_file(argv[n], flags);
  }
  return EXIT_SUCCESS;
}
