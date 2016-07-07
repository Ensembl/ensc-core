#include "sam.h"
#define MY_FUSEDFLAG 32768
/*
 print out a bam1_t entry, particularly the flags (for debugging)
 /nfs/ensembl/searle/progs/production_code/ensembl-trunk_0713/samtools_parallel
*/
void printBam(FILE *fp, bam1_t *b, bam_hdr_t *header) {    
  fprintf(fp, "%s %s %d %d %s (%d) %d %d %d\t\tP %d PP %d U %d MU %d R %d MR %d R1 %d R2 %d S %d QC %d D %d U %d\n",
                                  bam_get_qname(b), 
                                  header->target_name[b->core.tid], 
                                  b->core.pos, 
                                  //bam_calend(&b->core,bam1_cigar(b)),
                                  bam_endpos(b),
                                  header->target_name[b->core.mtid], 
                                  b->core.mtid, 
                                  b->core.mpos, 
                                  b->core.isize, 
                                  bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b)),
                                  b->core.flag & BAM_FPAIRED,
                                  b->core.flag & BAM_FPROPER_PAIR ? 1 : 0,
                                  b->core.flag & BAM_FUNMAP ? 1 : 0,
                                  b->core.flag & BAM_FMUNMAP ? 1 : 0,
                                  b->core.flag & BAM_FREVERSE ? 1 : 0,
                                  b->core.flag & BAM_FMREVERSE ? 1 : 0,
                                  b->core.flag & BAM_FREAD1 ? 1 : 0,
                                  b->core.flag & BAM_FREAD2 ? 1 : 0,
                                  b->core.flag & BAM_FSECONDARY ? 1 : 0,
                                  b->core.flag & BAM_FQCFAIL ? 1 : 0,
                                  b->core.flag & BAM_FDUP ? 1 : 0,
                                  b->core.flag & MY_FUSEDFLAG ? 1 : 0
                                  );
  fflush(fp);
}

long long countReadsInFile(char *inFName) {
  htsFile *in = hts_open(inFName, "rb");
  long long nUsableReads = 0;

  if (in == 0) {
    fprintf(stderr, "Fail to open BAM file %s\n", inFName);
    exit(1);
  }

  bam1_t *b = bam_init1();

  int cnt = 0;
  while (bam_read1(in->fp.bgzf, b) > 0 && b->core.tid >= 0) {
    if (!(b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP))) {
      nUsableReads++;
    }
    cnt++;
    if (!(cnt%1000000)) {
      fprintf(stderr,".");
      fflush(stdout);
    }
  }
  fprintf(stderr, "\n");

  hts_close(in);

  return nUsableReads;
}

