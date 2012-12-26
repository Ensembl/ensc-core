#include <stdio.h>

#include "EnsC.h"

#include "DBAdaptor.h"
#include "BaseAdaptor.h"
#include "Basic/Vector.h"
#include "Slice.h"
#include "StrUtil.h"



#include "sam.h"
#include "bam.h"

typedef struct {
  int beg, end;
  samfile_t *in;
} tmpstruct_t;

typedef struct mappingStruct {
  Slice *sourceSlice;
  Slice *destSlice;
  int    ori;
} Mapping;


Vector *getDestinationSlices(DBAdaptor *dba, char *assName);
Vector *getMappings(DBAdaptor *dba, Slice *destSlice, char *fromAssName, char *toAssName);
samfile_t *writeBamHeader(char *inFName, char *outFName, Vector *destinationSlices);
bam_header_t *bam_header_dup(const bam_header_t *h0);
int mapBam(char *fName, samfile_t *out, Mapping *mapping);
int mapLocation(Mapping *mapping, int pos);
void Bammap_usage();

int totOffEnd = 0;

int main(int argc, char *argv[]) {
  DBAdaptor *dba;
  StatementHandle *sth;
  ResultRow *row;
  Vector *destinationSlices;
  Vector *mappings;
  samfile_t *out;

  char *gtStr = ">";
  char *emptyStr = "";
  int argNum = 1;

  char *inFName  = NULL;
  char *outFName = NULL;

  char *dbHost = "ens-staging.internal.sanger.ac.uk";
  char *dbUser = "ensro";
  char *dbPass = NULL;
  char *dbName = "homo_sapiens_core_70_37";
  int   dbPort = 3306;

  char *sourceName = "GRCh37";
  char *destName   = "NCBI36";

  initEnsC();

//  primary_parser = secondary_parser = secondary_multi_parser = emptyStr;
// line_prefix = secondary_line_prefix = record_prefix = gtStr;

  while (argNum < argc) {
    char *arg = argv[argNum];
    char *val;

    if (argNum == argc-1) {
      Bammap_usage();
    }

    val = argv[++argNum];
//    printf("%s %s\n",arg,val);

    if (!strcmp(arg, "-i") || !strcmp(arg,"--in_file")) {
      StrUtil_copyString(&inFName,val,0);
    } else if (!strcmp(arg, "-o") || !strcmp(arg,"--out_file")) {
      StrUtil_copyString(&outFName,val,0);
    } else if (!strcmp(arg, "-h") || !strcmp(arg,"--host")) {
      StrUtil_copyString(&dbHost,val,0);
    } else if (!strcmp(arg, "-p") || !strcmp(arg,"--password")) {
      StrUtil_copyString(&dbPass,val,0);
    } else if (!strcmp(arg, "-P") || !strcmp(arg,"--port")) {
      dbPort = atoi(val);
    } else if (!strcmp(arg, "-n") || !strcmp(arg,"--name")) {
      StrUtil_copyString(&dbName,val,0);
    } else if (!strcmp(arg, "-s") || !strcmp(arg,"--source_ass")) {
      StrUtil_copyString(&sourceName,val,0);
    } else if (!strcmp(arg, "-d") || !strcmp(arg,"--dest_ass")) {
      StrUtil_copyString(&destName,val,0);
    } else {
      printf("Error in command line at %s\n\n",arg);
      Bammap_usage();
    }

    argNum++;
  }

  printf("Program for mapping BAM files between assemblies\n"
         "Steve M.J. Searle.  searle@sanger.ac.uk  Last update Dec 2012.\n");

  if (!inFName || !outFName) {
    Bammap_usage();
  }

  printf("Opening connection ...");
  dba = DBAdaptor_new(dbHost,dbUser,dbPass,dbName,dbPort,NULL);
  printf(" Done\n");

  destinationSlices = getDestinationSlices(dba, destName);

  out = writeBamHeader(inFName,outFName,destinationSlices);

  int i;
  for (i=0; i<Vector_getNumElement(destinationSlices); i++) {
    Slice *slice = Vector_getElementAt(destinationSlices,i);

    //if (!strcmp(Slice_getChrName(slice),"3")) break;

    mappings = getMappings(dba,slice,sourceName,destName);
    int j;
    for (j=0;j<Vector_getNumElement(mappings); j++) {
      Mapping *mapping = Vector_getElementAt(mappings,j);
      printf("%s\t%d\t%d\t%s\t%d\t%d\t%d\n",Slice_getChrName(mapping->destSlice),
                                            Slice_getChrStart(mapping->destSlice),
                                            Slice_getChrEnd(mapping->destSlice),
                                            Slice_getChrName(mapping->sourceSlice),
                                            Slice_getChrStart(mapping->sourceSlice),
                                            Slice_getChrEnd(mapping->sourceSlice),
                                            mapping->ori);
      
//      mapBam("ftp://ngs.sanger.ac.uk/scratch/project/searle/bams/GM12878_fixed_sorted.bam",out, mapping->sourceSlice);
      //if (Slice_getChrStart(mapping->destSlice)  > 97000000 && Slice_getChrStart(mapping->destSlice) < 98000000) {
      mapBam(inFName, out, mapping);
      //}
    }
  }

  printf("total off end = %d\n",totOffEnd);

  samclose(out);
  return 0;
}


void Bammap_usage() {
  printf("bammap \n"
         "  -i --in_file    Input BAM file to map from\n"
         "  -o --out_file   Output BAM file to write\n"
         "  -h --host       Database host name for db containing mapping\n"
         "  -n --name       Database name for db containing mapping\n"
         "  -p --password   Database password\n"
         "  -P --port       Database port\n"
         "  -s --source_ass Assembly name to map from\n"
         "  -d --dest_ass   Assembly name to map to\n");
  exit(1);
}

/* 
  Take header from input BAM file, replace the seq regions, and then write it out
*/
samfile_t *writeBamHeader(char *inFName, char *outFName, Vector *destinationSlices) {
  tmpstruct_t tmp;
  bam_header_t *destHeader;
  char out_mode[5];
  samfile_t *out;
  char line[1024];
  char *buff = NULL;

  char **tokens = NULL;
  int ntoken;

  tmp.in = samopen(inFName, "rb", 0);
  if (tmp.in == 0) {
    fprintf(stderr, "Fail to open BAM file %s\n", inFName);
    return NULL;
  }

  strcpy(out_mode, "wh");
  samopen("-",out_mode,tmp.in->header);

  /* Create the @SQ lines for the destination set of sequences */
  int i;
  for (i=0; i<Vector_getNumElement(destinationSlices); i++) {
    Slice *slice = Vector_getElementAt(destinationSlices,i);
    sprintf(line,"@SQ\tSN:%s\tLN:%d\n", Slice_getChrName(slice), 
                                        Slice_getChrEnd(slice));
    if (buff) {
      buff = StrUtil_appendString(buff,line);
    } else {
      buff = StrUtil_copyString(&buff,line,0);
    }
  }

  destHeader = bam_header_init();

  /* add non @SQ header lines from source header into destination header */
  StrUtil_tokenizeByDelim(&tokens, &ntoken, tmp.in->header->text, '\n');
  printf("ntoken = %d\n",ntoken);
  for (i=0;i<ntoken;i++) {
    printf("token = %s\n",tokens[i]);
    if (strncmp(tokens[i],"@SQ",3)) {
      /* HD must come first */
      if (!strncmp(tokens[i],"@HD",3)) {
        char *tmpBuff;
        StrUtil_copyString(&tmpBuff,tokens[i],0);
        tmpBuff = StrUtil_appendString(tmpBuff,"\n");
        tmpBuff = StrUtil_appendString(tmpBuff,buff);
        free(buff);
        buff = tmpBuff;
      } else {
        buff = StrUtil_appendString(buff,tokens[i]);
        buff = StrUtil_appendString(buff,"\n");
      }
    }
  }

  destHeader->l_text = strlen(buff);
  destHeader->text   = buff;
  printf("header->text = %s\n",buff);
  sam_header_parse(destHeader);

  strcpy(out_mode, "wbh");
  out = samopen(outFName,out_mode,destHeader);

  bam_init_header_hash(out->header);

  samclose(tmp.in);

  return out;
}


int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 9, 14, 1, 6, 5, 13, 3, 11, 7, 15 };

int mapBam(char *fName, samfile_t *out, Mapping *mapping) {
  tmpstruct_t tmp;
  int ref;
  char region[1024];
  int bCount;


  tmp.in = samopen(fName, "rb", 0);
  if (tmp.in == 0) {
    fprintf(stderr, "Fail to open BAM file %s\n", fName);
    return 1;
  }

  bam_index_t *idx;

  idx = bam_index_load(fName); // load BAM index
  if (idx == 0) {
    fprintf(stderr, "BAM indexing file is not available.\n");
    return 1;
  }
  sprintf(region,"%s:%d-%d", Slice_getChrName(mapping->sourceSlice), 
                             Slice_getChrStart(mapping->sourceSlice), 
                             Slice_getChrEnd(mapping->sourceSlice));
  bam_parse_region(tmp.in->header, region, &ref, &tmp.beg, &tmp.end);
  if (ref < 0) {
    fprintf(stderr, "Invalid region %s %d %d\n", Slice_getChrName(mapping->sourceSlice), 
                                                 Slice_getChrStart(mapping->sourceSlice), 
                                                 Slice_getChrEnd(mapping->sourceSlice));
    return 1;
  }

  bCount = 0;
  bam_iter_t iter = bam_iter_query(idx, ref, tmp.beg, tmp.end);
  bam1_t *b = bam_init1();
  int nOffEnds = 0;
  int nMateNonLocal = 0;

//  fprintf(stderr,"HERE slice chr name = %s\n", Slice_getChrName(mapping->destSlice));
//  samopen("-","wh",out->header);

  int32_t newtid = bam_get_tid(out->header, Slice_getChrName(mapping->destSlice));


  while (bam_iter_read(tmp.in->x.bam, iter, b) >= 0) {
    int end;

    bCount++;

    end = bam_calend(&b->core, bam1_cigar(b));
    if (end > tmp.end || b->core.pos < tmp.beg) {
      nOffEnds++;
      continue;
    }
    
    b->core.tid = newtid;

    if (b->core.mtid >= 0) {
      int32_t newmtid = bam_get_tid(out->header, tmp.in->header->target_name[b->core.mtid]);
      b->core.mtid = newmtid;

      if (!(b->core.mpos <= tmp.end && b->core.mpos >= tmp.beg && newmtid == newtid)) {
        nMateNonLocal++;
      } else {
        b->core.mpos = mapLocation(mapping, b->core.mpos+1) - 1;
      }
    }

/* pos should be left most position on reference so is either pos or end depending on mapping orientation */
    if (mapping->ori == 1) {
      b->core.pos = mapLocation(mapping, b->core.pos+1) - 1;
    } else {
      b->core.pos = mapLocation(mapping, end) - 1;
    }

    /* Switch strand if reverse orientation mapping */
    /* Also have to revcom the sequence */
    /* Also have to reverse the quality */
    /* Also have to reverse the cigar */
    if (mapping->ori == -1) {
      printf(" flag before = %d\n",b->core.flag);
      b->core.flag ^= BAM_FREVERSE;
      printf(" flag after = %d\n",b->core.flag);

      uint8_t *seq;
      seq = bam1_seq(b);
      int i;
      uint8_t *revseq = calloc(1,(b->core.l_qseq+1)/2);
      
/* print out existing sequence */
/*
      printf("forward: ");
      for (i=0;i<b->core.l_qseq;i++) {
        char base = bam_nt16_rev_table[bam1_seqi(seq, i)];
        printf("%c", base);
      }
      printf("\n");
*/
/* make rev com seq */
      for (i=0;i<b->core.l_qseq;i++) {
/*
        char base = bam_nt16_rev_table[bam1_seqi(seq, b->core.l_qseq-i-1)];
        revseq[i/2] |= bam_nt16_table[(int)base] << 4*(1-i%2);
*/
        uint8_t ib = seq_comp_table[bam1_seqi(seq, b->core.l_qseq-i-1)];
        revseq[i/2] |= ib << 4*(1-i%2);
      }
/* print rev com seq */
/*
      printf("reverse: ");
      for (i=0;i<b->core.l_qseq;i++) {
        char base = bam_nt16_rev_table[bam1_seqi(revseq, i)];
        printf("%c", base);
      }
      printf("\n");
*/
      memcpy(seq,revseq,(b->core.l_qseq+1)/2);
      free(revseq);

      uint8_t *qual;
      qual = bam1_qual(b);
      uint8_t *revqual = calloc(1,b->core.l_qseq);
      for (i=0;i<b->core.l_qseq;i++) {
        revqual[i] = qual[b->core.l_qseq-i-1];
      }
      memcpy(qual,revqual,b->core.l_qseq);
      free(revqual);
      
      int *revcig = calloc(1,b->core.n_cigar*sizeof(int));
      int *cig = bam1_cigar(b);
      
      for (i=0;i<b->core.n_cigar;i++) {
        revcig[i] = cig[b->core.n_cigar-1-i];
      }
      memcpy(cig,revcig,b->core.n_cigar * 4);

      free(revcig);
    }

/* recalculate bin, not sure if need this for all coords, or just - ori ones */
    end = bam_calend(&b->core, bam1_cigar(b));
    b->core.bin = bam_reg2bin(b->core.pos,end);

    if (!bam_write1(out->x.bam, b)) {
      printf("Failed writing bam entry\n");
    }
  }
  printf("Number of reads (bam_iter)  = %d  num off ends %d num non local mates %d\n",bCount,nOffEnds,nMateNonLocal);
  totOffEnd+=nOffEnds;

  bam_iter_destroy(iter);
  bam_destroy1(b);

  bam_index_destroy(idx);

  samclose(tmp.in);
  return 0;
}


int mapLocation(Mapping *mapping, int pos) {
  int fromStart  = Slice_getChrStart(mapping->sourceSlice);
  int fromEnd    = Slice_getChrEnd(mapping->sourceSlice);
  int toStart    = Slice_getChrStart(mapping->destSlice);
  int toEnd      = Slice_getChrEnd(mapping->destSlice);

  if (pos < fromStart || pos > fromEnd) {
    printf("Error: tried to map position out of range pos = %d range = %s %d %d\n", pos, Slice_getChrName(mapping->sourceSlice), fromStart, fromEnd);
  }

  if (mapping->ori == 1) {
    pos += (toStart-fromStart);
  } else {

    int offset = pos-fromStart;
    //printf("pos before = %d\n",pos);
    //printf("offset = %d\n",offset);
    pos = toEnd-offset;
    //printf("pos after  = %d\n\n",pos);
  }

  return pos;
}


Vector *getMappings(DBAdaptor *dba, Slice *destSlice, char *fromAssName, char *toAssName) {
  StatementHandle *sth;
  ResultRow *row;
  char qStr[1024];
  char *coordSys1;
  char *coordSys2;
  
  Vector *mappingVector = Vector_new();

  int i;
  for (i=0;i<2;i++) {
    if (!i) {
      coordSys1 = toAssName;
      coordSys2 = fromAssName;
    } else {
      coordSys1 = fromAssName;
      coordSys2 = toAssName;
    }
    sprintf(qStr,
            "select sr1.name, a.asm_start, a.asm_end, sr2.name, a.cmp_start, a.cmp_end, a.ori"
            " from seq_region sr1, seq_region sr2, assembly a, coord_system cs1, coord_system cs2"
            " where sr1.seq_region_id=a.asm_seq_region_id and"
            "       sr2.seq_region_id=a.cmp_seq_region_id and"
            "       cs1.coord_system_id=sr1.coord_system_id and"
            "       cs2.coord_system_id=sr2.coord_system_id and"
            "       cs1.version='%s' and"
            "       cs2.version='%s' and"
            "       sr%d.name='%s' order by a.%s_start",
            coordSys1, coordSys2, i+1,  Slice_getChrName(destSlice), (!i ? "asm" : "cmp"));
  
    printf("%s\n",qStr);
  
    sth = dba->dbc->prepare(dba->dbc,qStr,strlen(qStr));
  
    sth->execute(sth);
  
  
    while (row = sth->fetchRow(sth)) {
      char *destName;
      int   destStart;
      int   destEnd;
      char *sourceName;
      int   sourceStart;
      int   sourceEnd;
      int   ori;
      Mapping *mapping;
      Slice *dest;
      Slice *source;
  
      if (row->col(row,0))  destName   = row->getStringAt(row,0);
      if (row->col(row,1))  destStart  = row->getIntAt(row,1);
      if (row->col(row,2))  destEnd    = row->getIntAt(row,2);
      if (row->col(row,3))  sourceName   = row->getStringAt(row,3);
      if (row->col(row,4))  sourceStart  = row->getIntAt(row,4);
      if (row->col(row,5))  sourceEnd    = row->getIntAt(row,5);
      if (row->col(row,6))  ori          = row->getIntAt(row,6);
       
      dest    = Slice_new(destName,destStart,destEnd,1,toAssName,NULL,0,0);
      source  = Slice_new(sourceName,sourceStart,sourceEnd,1,fromAssName,NULL,0,0);
      mapping = calloc(1,sizeof(Mapping)); 
  
      if (!i) {
        mapping->sourceSlice = source;
        mapping->destSlice   = dest;
        mapping->ori         = ori;
      } else {
        mapping->sourceSlice = dest;
        mapping->destSlice   = source;
        mapping->ori         = ori;
      }
      
      Vector_addElement(mappingVector,mapping); 
  
    }
  
    sth->finish(sth);
  }

  return mappingVector;
}

Vector *getDestinationSlices(DBAdaptor *dba, char *assName) {
  StatementHandle *sth;
  ResultRow *row;
  char qStr[1024];
  char *assemblyType = DBAdaptor_getAssemblyType(dba);

  if (!strcmp(assemblyType,assName)) {
    sprintf(qStr,
            "select distinct(sr.name), sr.length, sr.coord_system_id"
            " from seq_region sr,"
            "      seq_region_attrib sra,"
            "      attrib_type at,"
            "      coord_system cs "
            " where cs.version = '%s' and"
            "       cs.coord_system_id=sr.coord_system_id and"
            "       sra.seq_region_id=sr.seq_region_id and"
            "       sra.attrib_type_id=at.attrib_type_id and"
            "       at.code='toplevel' and"
            "       sr.seq_region_id not in "
            "            (select sr1.seq_region_id from seq_region sr1, seq_region_attrib sra1, attrib_type at1"
            "             where sr1.seq_region_id=sra1.seq_region_id and"
            "                   sra1.attrib_type_id=at1.attrib_type_id and"
            "                   at1.code = 'non_ref') order by length desc",
            assName);
  } else {
    sprintf(qStr,
            "select distinct(sr.name), sr.length, sr.coord_system_id"
            " from seq_region sr,"
            "      coord_system cs "
            " where cs.version = '%s' and"
            "       cs.coord_system_id=sr.coord_system_id order by sr.length desc",
            assName);
  }

  sth = dba->dbc->prepare(dba->dbc,qStr,strlen(qStr));

  sth->execute(sth);

  Vector *toplevelSliceVector = Vector_new();

  while (row = sth->fetchRow(sth)) {
    char *name;
    int   length;
    int   cs;
    Slice *slice;

    if (row->col(row,0))  name   = row->getStringAt(row,0);
    if (row->col(row,1))  length = row->getIntAt(row,1);
    if (row->col(row,2))  cs     = row->getIntAt(row,2);
     
    slice = Slice_new(name,1,length,1,assName,NULL,0,0);

    Vector_addElement(toplevelSliceVector,slice); 

  }


  sth->finish(sth);

  return toplevelSliceVector;
}
