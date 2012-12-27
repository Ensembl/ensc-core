#include <stdio.h>

#include "EnsC.h"

#include "DBAdaptor.h"
#include "BaseAdaptor.h"
#include "Basic/Vector.h"
#include "Slice.h"
#include "StrUtil.h"

#include "sam.h"
#include "bam.h"

typedef struct mappingStruct {
  Slice *sourceSlice;
  Slice *destSlice;
  int    ori;
} Mapping;

typedef struct readMapStatsStruct {
  int nRead;
  int nWritten;
  int nOverEnds;
  int nRemoteMate;
  int nReversed;
  int nUnmappedMate;
} ReadMapStats;

Vector *      getDestinationSlices(DBAdaptor *dba, char *assName);
Vector *      getMappings(DBAdaptor *dba, char *seqName, char *fromAssName, char *toAssName, int rev);
samfile_t *   writeBamHeader(char *inFName, char *outFName, Vector *destinationSlices);
bam_header_t *bam_header_dup(const bam_header_t *h0);
int           mapBam(char *fName, samfile_t *out, Mapping *mapping, ReadMapStats *regionStats, samfile_t *in, bam_index_t *idx, Vector **mappingVectors);
int           mapLocation(Mapping *mapping, int pos);
void          Bammap_usage();
int           mapRemoteLocation(Vector **mappingVectors, int seqid, int pos);
Vector **     getMappingVectorsBySourceRegion(DBAdaptor *dba, samfile_t *in, char *sourceName, char *destName);

int main(int argc, char *argv[]) {
  DBAdaptor *dba;
  StatementHandle *sth;
  ResultRow *row;
  Vector *destinationSlices;
  Vector *mappings;
  samfile_t *out;

  ReadMapStats totalStats;
  ReadMapStats regionStats;

  memset(&totalStats, 0, sizeof(ReadMapStats));

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
    } else if (!strcmp(arg, "-u") || !strcmp(arg,"--user")) {
      StrUtil_copyString(&dbUser,val,0);
    } else if (!strcmp(arg, "-s") || !strcmp(arg,"--source_ass")) {
      StrUtil_copyString(&sourceName,val,0);
    } else if (!strcmp(arg, "-d") || !strcmp(arg,"--dest_ass")) {
      StrUtil_copyString(&destName,val,0);
    } else {
      fprintf(stderr,"Error in command line at %s\n\n",arg);
      Bammap_usage();
    }

    argNum++;
  }

  printf("Program for mapping BAM files between assemblies\n"
         "Steve M.J. Searle.  searle@sanger.ac.uk  Last update Dec 2012.\n");

  if (!inFName || !outFName) {
    Bammap_usage();
  }

  //printf("Opening connection ...");
  dba = DBAdaptor_new(dbHost,dbUser,dbPass,dbName,dbPort,NULL);
  //printf(" Done\n");

  destinationSlices = getDestinationSlices(dba, destName);

  out = writeBamHeader(inFName,outFName,destinationSlices);

  samfile_t *in = samopen(inFName, "rb", 0);
  if (in == 0) {
    fprintf(stderr, "Fail to open BAM file %s\n", inFName);
    return 1;
  }

  bam_init_header_hash(in->header);

  bam_index_t *idx;
  idx = bam_index_load(inFName); // load BAM index
  if (idx == 0) {
    fprintf(stderr, "BAM indexing file is not available.\n");
    return 1;
  }

  // Load all mappings into an array of vectors arranged by source region id
  // Only used for remote mate mapping, where we don't know in advanced which
  // region the mapping will be to
  Vector **mappingVectors = getMappingVectorsBySourceRegion(dba, in, sourceName, destName);

  int i;
  for (i=0; i<Vector_getNumElement(destinationSlices); i++) {
    Slice *slice = Vector_getElementAt(destinationSlices,i);

    printf("Working on '%s'\n",Slice_getChrName(slice));
    //if (!strcmp(Slice_getChrName(slice),"3")) break;

    mappings = getMappings(dba,Slice_getChrName(slice),sourceName,destName,0);
    int j;
    for (j=0;j<Vector_getNumElement(mappings); j++) {
      Mapping *mapping = Vector_getElementAt(mappings,j);
//      printf("%s\t%d\t%d\t%s\t%d\t%d\t%d\n",Slice_getChrName(mapping->destSlice),
//                                            Slice_getChrStart(mapping->destSlice),
//                                            Slice_getChrEnd(mapping->destSlice),
//                                            Slice_getChrName(mapping->sourceSlice),
//                                            Slice_getChrStart(mapping->sourceSlice),
//                                            Slice_getChrEnd(mapping->sourceSlice),
//                                            mapping->ori);
      
//      mapBam("ftp://ngs.sanger.ac.uk/scratch/project/searle/bams/GM12878_fixed_sorted.bam",out, mapping->sourceSlice);
      //if (Slice_getChrStart(mapping->destSlice)  > 97000000 && Slice_getChrStart(mapping->destSlice) < 98000000) {
      memset(&regionStats, 0, sizeof(ReadMapStats));

      mapBam(inFName, out, mapping, &regionStats, in, idx, mappingVectors);

      totalStats.nRead         += regionStats.nRead;
      totalStats.nWritten      += regionStats.nWritten;
      totalStats.nReversed     += regionStats.nReversed;
      totalStats.nOverEnds     += regionStats.nOverEnds;
      totalStats.nRemoteMate   += regionStats.nRemoteMate;
      totalStats.nUnmappedMate += regionStats.nUnmappedMate;
      //}
    }
  }

  printf(" Total reads read in mapped regions           %d\n", totalStats.nRead);
  printf(" Total reads written in mapped regions        %d\n", totalStats.nWritten);
  printf(" Total reads where orientation reversed       %d\n", totalStats.nReversed);
  printf(" Total reads with remotely located mates      %d\n", totalStats.nRemoteMate);
  printf(" Total reads extending beyond ends of region  %d (not written)\n", totalStats.nOverEnds);
  printf(" Total reads with unmapped mates              %d (not written)\n", totalStats.nUnmappedMate);

  samclose(out);
  bam_index_destroy(idx);

  samclose(in);
  return 0;
}


void Bammap_usage() {
  printf("bammap \n"
         "  -i --in_file    Input BAM file to map from\n"
         "  -o --out_file   Output BAM file to write\n"
         "  -h --host       Database host name for db containing mapping\n"
         "  -n --name       Database name for db containing mapping\n"
         "  -u --user       Database user\n"
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
  bam_header_t *destHeader;
  char          out_mode[5];
  samfile_t *   in;
  samfile_t *   out;
  char          line[1024];
  char *        buff = NULL;
  char **       tokens = NULL;
  int           ntoken;

  in = samopen(inFName, "rb", 0);
  if (in == 0) {
    fprintf(stderr, "Fail to open BAM file %s\n", inFName);
    return NULL;
  }

  //strcpy(out_mode, "wh");
  //samopen("-",out_mode,in->header);

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
  StrUtil_tokenizeByDelim(&tokens, &ntoken, in->header->text, '\n');
  //printf("ntoken = %d\n",ntoken);
  for (i=0;i<ntoken;i++) {
    //printf("token = %s\n",tokens[i]);
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
  //printf("header->text = %s\n",buff);
  sam_header_parse(destHeader);

  strcpy(out_mode, "wbh");
  out = samopen(outFName,out_mode,destHeader);

/* Need this hash initialised for looking up tids */
  bam_init_header_hash(out->header);

  samclose(in);

  return out;
}

int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 9, 14, 1, 6, 5, 13, 3, 11, 7, 15 };

int mapBam(char *fName, samfile_t *out, Mapping *mapping, ReadMapStats *regionStats, samfile_t *in, bam_index_t *idx, Vector **mappingVectors) {
  int  ref;
  int  begRange;
  int  endRange;
  char region[1024];

  sprintf(region,"%s:%d-%d", Slice_getChrName(mapping->sourceSlice), 
                             Slice_getChrStart(mapping->sourceSlice), 
                             Slice_getChrEnd(mapping->sourceSlice));
  bam_parse_region(in->header, region, &ref, &begRange, &endRange);
  if (ref < 0) {
    fprintf(stderr, "Invalid region %s %d %d\n", Slice_getChrName(mapping->sourceSlice), 
                                                 Slice_getChrStart(mapping->sourceSlice), 
                                                 Slice_getChrEnd(mapping->sourceSlice));
    return 1;
  }

  bam_iter_t iter = bam_iter_query(idx, ref, begRange, endRange);
  bam1_t *b = bam_init1();

//  fprintf(stderr,"HERE slice chr name = %s\n", Slice_getChrName(mapping->destSlice));
//  samopen("-","wh",out->header);

  int32_t newtid = bam_get_tid(out->header, Slice_getChrName(mapping->destSlice));

  while (bam_iter_read(in->x.bam, iter, b) >= 0) {
    int end;

    regionStats->nRead++;

    end = bam_calend(&b->core, bam1_cigar(b));

    if (end > endRange || b->core.pos < begRange) {
      regionStats->nOverEnds++;
      continue;
    }
    
    b->core.tid = newtid;

    if (b->core.mtid >= 0) {
      int32_t newmtid = bam_get_tid(out->header, in->header->target_name[b->core.mtid]);

      if (!(b->core.mpos <= endRange && b->core.mpos >= begRange && newmtid == newtid)) {
        regionStats->nRemoteMate++;
        if ((b->core.mpos = mapRemoteLocation(mappingVectors, b->core.mtid, b->core.mpos+1) - 1) < 0) {
          regionStats->nUnmappedMate++;
          continue;
        }
      } else {
        if ((b->core.mpos = mapLocation(mapping, b->core.mpos+1) - 1) < 0) {
          regionStats->nUnmappedMate++;
          continue;
        }
      }

      b->core.mtid = newmtid;
    }

/* pos should be left most position on reference so is either pos or end depending on mapping orientation */
    if (mapping->ori == 1) {
      b->core.pos = mapLocation(mapping, b->core.pos+1) - 1;
    } else {
      b->core.pos = mapLocation(mapping, end) - 1;
    }

    /* toggle (XOR) rev com flag if reverse orientation mapping */
    /* Also have to revcom the sequence */
    /* Also have to reverse the quality */
    /* Also have to reverse the cigar */
    if (mapping->ori == -1) {
      b->core.flag ^= BAM_FREVERSE;

      regionStats->nReversed++;

      uint8_t *seq;
      uint8_t *qual;
      seq = bam1_seq(b);
      qual = bam1_qual(b);

      int i;
      uint8_t *revseq  = calloc(1, (b->core.l_qseq+1)/2 + b->core.l_qseq);
      uint8_t *revqual = revseq + (b->core.l_qseq+1)/2;
      
/* make rev com seq and reverse qual at same time*/
      for (i=0;i<b->core.l_qseq;i++) {
        uint8_t ib = seq_comp_table[bam1_seqi(seq, b->core.l_qseq-i-1)];
        revseq[i/2] |= ib << 4*(1-i%2);
        revqual[i] = qual[b->core.l_qseq-i-1];
      }

      memcpy(seq, revseq, (b->core.l_qseq+1)/2 + b->core.l_qseq);
      free(revseq);

/*
      uint8_t *revqual = calloc(1,b->core.l_qseq);
      for (i=0;i<b->core.l_qseq;i++) {
        revqual[i] = qual[b->core.l_qseq-i-1];
      }
      memcpy(qual,revqual,b->core.l_qseq);
      free(revqual);
*/
      
      int *revcig = calloc(1,b->core.n_cigar*4);
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

    regionStats->nWritten++;

    if (!bam_write1(out->x.bam, b)) {
      fprintf(stderr, "Failed writing bam entry\n");
    }
  }
  //printf("Number of reads read = %d  num off ends %d num non local mates %d\n",regionStats->nRead,
  //                                                                             regionStats->nOverEnds,
  //                                                                             regionStats->nRemoteMate);

  bam_iter_destroy(iter);
  bam_destroy1(b);

  return 0;
}

// Not a very efficient routine
// Could do a binary search into the vector to speed up finding the right mapping
//    Not many reads needed this call so I didn't bother
int mapRemoteLocation(Vector **mappingVectors, int seqid, int pos) {
  Mapping *mapping = NULL;
  Vector  *mapVec = mappingVectors[seqid];
  int i;

  /* Find the mapping containing the location */
  for (i=0; i<Vector_getNumElement(mapVec); i++) {
    Mapping *m = Vector_getElementAt(mapVec,i);
    char *fromChr = Slice_getChrName(m->sourceSlice);
    int fromStart = Slice_getChrStart(m->sourceSlice);
    int fromEnd   = Slice_getChrEnd(m->sourceSlice);

    //printf("Comparing %d to slice %s %d %d\n",pos,fromChr,fromStart,fromEnd);

    if (pos >= fromStart && pos <= fromEnd) {
      mapping = m;
      //printf("Match %d to slice %s %d %d\n",pos,fromChr,fromStart,fromEnd);
      break;
    }
  }

  if (!mapping) {
    //fprintf(stderr,"Mate location lies outside any mapped region pos = %d\n", pos);
    return -1;
  }

  /* call mapLocation */
  return mapLocation(mapping, pos);
}

int mapLocation(Mapping *mapping, int pos) {
  int fromStart  = Slice_getChrStart(mapping->sourceSlice);
  int fromEnd    = Slice_getChrEnd(mapping->sourceSlice);
  int toStart    = Slice_getChrStart(mapping->destSlice);
  int toEnd      = Slice_getChrEnd(mapping->destSlice);

  if (pos < fromStart || pos > fromEnd) {
    fprintf(stderr,"Error: tried to map position out of range pos = %d range = %s %d %d\n", pos, Slice_getChrName(mapping->sourceSlice), fromStart, fromEnd);
    return -1;
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


Vector **getMappingVectorsBySourceRegion(DBAdaptor *dba, samfile_t *in, char *sourceName, char *destName) {
  int i;

  Vector **mappingVectors = calloc(in->header->n_targets, sizeof(Vector *));

  for (i=0;i<in->header->n_targets;i++) {
    char *seqName = in->header->target_name[i];

    // Note reverse mapping direction to key on source
    // Use flag to getMappings to fill slices with correct direction
    mappingVectors[i] = getMappings(dba,seqName,destName,sourceName, 1); 
  }

  return mappingVectors;
}

// rev flag is a bit of a hack, to enable fetching by source region, but switching slices so
// it looks like the mapping is the other way round 
Vector *getMappings(DBAdaptor *dba, char *seqName, char *fromAssName, char *toAssName, int rev) {
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
            coordSys1, coordSys2, i+1,  seqName, (!i ? "asm" : "cmp"));
  
    //printf("%s\n",qStr);
  
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
  
      if (!rev) {
        if (!i) {
          mapping->sourceSlice = source;
          mapping->destSlice   = dest;
          mapping->ori         = ori;
        } else {
          mapping->sourceSlice = dest;
          mapping->destSlice   = source;
          mapping->ori         = ori;
        }
      } else {
        if (!i) {
          mapping->sourceSlice = dest;
          mapping->destSlice   = source;
          mapping->ori         = ori;
        } else {
          mapping->sourceSlice = source;
          mapping->destSlice   = dest;
          mapping->ori         = ori;
        }
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
