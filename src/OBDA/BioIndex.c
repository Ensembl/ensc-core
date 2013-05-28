/****************************************************************\
*                                                                *
*  Library for finding sequences using ODBA format indexes       *
*                                                                *
*  Steve Searle         mailto:searle@sanger.ac.uk               *
*  Copyright (C) 2000,2001.  All Rights Reserved.                *
*                                                                *
*  This source code is distributed under the terms of the        *
*  GNU Lesser General Public License. See the file COPYING       *
*  or http://www.fsf.org/copyleft/lesser.html for details        *
*                                                                *
*  If you use this code, please keep this notice intact.         *
*                                                                *
\****************************************************************/

#include "BioIndex.h"
#include "Vector.h"
#include "StrUtil.h"
#include "Error.h"
#include "EnsC.h"
#include "Sequence.h"

#include <ctype.h> /* for isspace() */
#include <string.h>

/* for scandir etc. */
#include <dirent.h>
/* for stat */
#include <sys/stat.h>
/* for MAXPATHLEN */
#include <sys/param.h>
 
int selectPrefix(const struct dirent *DirEnt);
char *currentPrefix;  
#define MAXSTRLEN 1024


int areintquiet(const char *string);

/* ---------------------------------------------------------------- */

static void getFullFName(char *full, char *path, char *fname) {
  sprintf(full,"%s/%s",path,fname);
} 


static void *xcalloc(int size, int nelem) {
  void *ptr;

  if ((ptr = calloc(size,nelem)) == NULL) {
    
    printf("Error:Failed in xcalloc\n");
    exit(1);
  }
  return ptr;
}

static void *xrealloc(void *ptr, int size) {
  void *tmp_ptr;

  if ((tmp_ptr = realloc(ptr,size)) == NULL) {

    printf("Error:Failed in xrealloc\n");
    exit(1);
  }
  return tmp_ptr;
}


static BioIndex_IndexFile *BioIndex_IndexFile_create(char *fname) {
  BioIndex_IndexFile *bif = xcalloc(sizeof(BioIndex_IndexFile),1);
  char header[1024];
  int ntoread;
  int i;
                      
  if ((bif->handle = fopen(fname,"r")) == NULL) {
    Error_write(EOBDA,"BioIndex_IndexFile_create",ERR_SEVERE,"Failed opening index file %s\n", fname);
  }

  /* Read header. "%4.4d" (record width)  + nfield x " %4.4d" */
  if (fread(header,1,4,bif->handle) != 4) { 
    Error_write(EOBDA,"BioIndex_IndexFile_create",ERR_SEVERE,"Failed reading header from index file %s\n", fname);
  }

  sscanf(header,"%d",&(bif->rec_length));
  printf("Rec length = %d\n",bif->rec_length);
/* TODO Add checks */
    
  bif->header_len = ftell(bif->handle);

  fseek(bif->handle, 0, SEEK_END);

  bif->num_record = (ftell(bif->handle) - bif->header_len) / bif->rec_length;          
  return bif; 
}

static void BioIndex_IndexFile_destroy(BioIndex_IndexFile *bif) {
  Error_write(EOBDA,"BioIndex_IndexFile_destroy",ERR_SEVERE,"Not implemented");
}


static void BioIndex_open_indices(BioIndex *bi, char *path) {
  char *primary_fname = xcalloc(MAXPATHLEN+1,1);
  char  line[1024];
  int   i;
  char  FName[MAXPATHLEN];
  char  FullFName[MAXPATHLEN];

  sprintf(primary_fname,"%s/key_%s.key",path,bi->primary_namespace);
  bi->primary_index = BioIndex_IndexFile_create(primary_fname);
   
  bi->secondary_indices = xcalloc(sizeof(BioIndex_IndexFile *),bi->num_secondary);

  for (i=0; i<bi->num_secondary; i++) {
    sprintf(FName,"id_%s.index",bi->secondary_namespaces[i]);
    printf("Secondary = %s\n",FName);
    getFullFName(FullFName, path, FName);
    bi->secondary_indices[i] = BioIndex_IndexFile_create(FullFName);
  }   

  free(primary_fname);
}

static void BioIndex_addLocation(BioIndex *bi, BioIndex_Location **locationsP,
                        char *header, long start, 
                        long end, int index, BioIndex_Parser_Func parser) {
  char *id;
  BioIndex_Location *location;

  bi->primary_index->num_record++;
  if (!(bi->primary_index->num_record % 1000000)) {
    *locationsP = xrealloc((*locationsP), 
                           sizeof(BioIndex_Location) *(bi->primary_index->num_record + 1000000));
  }
  location = &(*locationsP)[bi->primary_index->num_record-1];

  if ((id = parser(header)) != NULL) {
    //printf("id = %s\n",id);
    //fflush(stdout);
    location->fileIndex = index;
    location->start = start;
    location->length= end-start+1;
    location->entry = id;
  }
}

static void BioIndex_close_indices(BioIndex *bi) {
  int i;

  fclose(bi->primary_index->handle);
  for (i=0; i<bi->num_secondary; i++) {
    fclose(bi->secondary_indices[i]->handle);
  } 
}

int BioIndex_get_key_from_line(char *line, char *entry) {
  char *LineP  = line;
  char *EntryP = entry;

  while (*LineP && *LineP != '\t') {
    *EntryP = *LineP;
    EntryP++;
    LineP++;
  }
  *EntryP = '\0';
  return 1;
}

int BioIndex_parse_primary_record(BioIndex *bi, 
                                  char *line, 
                                  BioIndex_Location *loc) {
  
  int  ntok;
  char id[1024];
  char tmp_string[1024];
  char *ch_p = line;  

/* 
  if (strlen(line) != bi->primary_index->rec_length) {
    g_error("ERROR: Failed parsing line %s (incorrect length)\n", line);
  }
*/
    
  if ((ntok = sscanf(line,"%hd %ld %d",&(loc->fileIndex), 
                     &(loc->start),&(loc->length))) != 3) {
    printf("N token = %d\n",ntok);
    Error_write(EOBDA,"BioIndex_parse_primary_record",ERR_SEVERE,"ERROR: Failed parsing primary record %s\n", line);
  } 
  return 1;
}

char * BioIndex_get_entry(BioIndex *bi, BioIndex_Location *loc) {
  FILE *FpDat;
  char  faFile[1024];
  BioIndex_FileID *fid =  Vector_getElementAt(bi->fileids,loc->fileIndex);
  char *data = xcalloc(loc->length+1,1);

  sprintf(faFile,"%s", fid->name);
  
  if ((FpDat = fopen(faFile,"r")) == NULL) {
    Error_write(EOBDA,"BioIndex_get_entry",ERR_SEVERE,"ERROR: Failed opening %s for read\n",faFile);
  }
  
  fseek(FpDat,loc->start,SEEK_SET);

  fread(data,1,loc->length,FpDat); 

/*
    printf("Data = %s",data);
*/
  
  fclose(FpDat);
  return data;
}

Sequence * BioIndex_get_fasta_sequence(BioIndex *bi, BioIndex_Location *loc,
                                       char *id) {
  char *data;
  Sequence *seq;
  char *ChP;
  char *RepP;

  data = BioIndex_get_entry(bi, loc);

  RepP = data;
  ChP = strchr(data,'\n');
  if (!ChP) {
    Error_write(EOBDA,"BioIndex_get_fasta_sequence",ERR_SEVERE,"Fasta format error for %s\n",loc->entry);
  }

  while (*ChP != '\0') {
    if (*ChP!='\n') {
      *RepP = *ChP;
      RepP++;
    }
    ChP++;
  }
  *RepP = '\0';

  seq = Sequence_new();
  Sequence_setName(seq,id);
  Sequence_setSeq(seq,data);
    
  free(data);

  return seq;
}

int BioIndex_get_fasta_entry(BioIndex *bi, BioIndex_Location *loc) {
  FILE *FpDat;
  char  faFile[1024];
  char  line[1024];
  BioIndex_FileID *fid =  Vector_getElementAt(bi->fileids,loc->fileIndex);

  sprintf(faFile,"%s", fid->name);
  
  if ((FpDat = fopen(faFile,"r")) == NULL) {
    Error_write(EOBDA,"BioIndex_get_fasta_entry",ERR_SEVERE,"ERROR: Failed opening %s for read\n",faFile);
  }
  
  fseek(FpDat,loc->start,SEEK_SET);
  if (fgets(line,1024,FpDat) == NULL) {
    Error_write(EFGETSNULL,"BioIndex_get_fasta_entry",ERR_SEVERE,"Failed reading first line of fasta entry %s from %s",loc->entry, faFile);
  }
  
  if (line[0] != '>' || !strstr(line, loc->entry)) {
    Error_write(EOBDA,"BioIndex_get_fasta_entry",ERR_SEVERE,"ERROR: Seek didn't go to correct entry for %s\n",loc->entry);
  }
  
  do {
    printf("%s",line);
    // Keep compiler happy by getting return value even though not really needed because of feof check just after
    char *ret = fgets(line,1024,FpDat);
  } while(!feof(FpDat) && line[0] != '>');
  
  fclose(FpDat);
  return 1;
}

static int BioIndex_count_strings(char **strings) {
  int count = 0;

  while (strings[count] != NULL) {
    count++;
  }
  return count;
}

static int BioIndex_binary_search_file(FILE *fp, int reclen, char *key,
                                       int nrec, int headlen, 
                                       long *recNum, 
                                       char *result_string){
  long   curRec  = nrec/2;
  long   lowRec  = 0;
  long   highRec = nrec;

  int   cmpVal;

  char  line[1024];
  char  curDBKey[1024];

  while (lowRec <= highRec) {
    fseek(fp, curRec*reclen+headlen, SEEK_SET);
    fread(line, reclen, 1, fp);
    line[reclen] = '\0';
/*
        printf("Looking at line %s\n",line);
*/
  
    BioIndex_get_key_from_line(line, curDBKey);
    cmpVal = strcmp(curDBKey, key);
/*
        printf("Compared |%s| to |%s|, got %d\n",curDBKey,key,cmpVal);
        printf("lowRec = %d, highRec = %d, curRec = %d\n",lowRec,highRec,curRec);
*/
    
    if (cmpVal == 0) {
      strcpy(result_string,line);
      *recNum = curRec;
      return TRUE;  
      break;
    } else if (cmpVal < 0) {
      lowRec = curRec+1;
      curRec = (highRec-lowRec)/2 + lowRec;
    } else {
      highRec = curRec-1;
      curRec = (highRec-lowRec)/2 + lowRec;
    }
  }
  if (cmpVal != 0) {
    Error_write(EOBDA,"BioIndex_binary_search_file",ERR_SEVERE,"ERROR: Failed finding |%s|\n", key);
  }
  return FALSE;
}



BioIndex_Location * BioIndex_get_by_primary_key(BioIndex *bi, char *key) {
  char * result_string = xcalloc(bi->primary_index->rec_length+1,1);
  BioIndex_Location * loc = NULL;
  long   recNum;

  if (BioIndex_binary_search_file(bi->primary_index->handle,
                                  bi->primary_index->rec_length, 
                                  key,
                                  bi->primary_index->num_record,
                                  bi->primary_index->header_len,
                                  &recNum,
                                  result_string)) {
    loc = xcalloc(sizeof(BioIndex_Location), 1);
    BioIndex_parse_primary_record(bi, &result_string[strlen(key)], loc);
  }

  free(result_string);
  return loc;
}

static char *BioIndex_get_secondary_value(char *result_string) {
  char * val_ptr;

  val_ptr = strchr(result_string,'\t');   
    
  if (val_ptr == NULL) {
    Error_write(EOBDA,"BioIndex_get_secondary_value",ERR_SEVERE,"Failed getting value from secondary index record %s\n",
            result_string);
  }

  return val_ptr+1;
}

void BioIndex_get_match_bounds(BioIndex *bi, char *key,
                               BioIndex_IndexFile *index,
                               long recNum,
                               long *start, long *end) {
  long   curRec;
  long   firstRec;
  char * result_string = xcalloc(index->rec_length+1,1);
  char   line[1024];
  char   curDBKey[1024];

  firstRec = curRec = recNum;
     
  while (curRec >= 0) {
    fseek(index->handle, curRec*index->rec_length+index->header_len, SEEK_SET);
    fread(line, index->rec_length, 1, index->handle);
    line[index->rec_length] = '\0';
  
    BioIndex_get_key_from_line(line, curDBKey);
    if (strcmp(curDBKey,key)) {
      break;
    }
    curRec--;
  }
  *start = curRec; 

  curRec = firstRec;
  while (curRec < index->num_record) {
    fseek(index->handle, curRec*index->rec_length+index->header_len, SEEK_SET);
    fread(line, index->rec_length, 1, index->handle);
    line[index->rec_length] = '\0';
  
    BioIndex_get_key_from_line(line, curDBKey);
    if (strcmp(curDBKey,key)) {
      break;
    }
    curRec++;
  }
  *end = curRec; 

  free(result_string);
}

Vector* BioIndex_get_by_secondary_key(BioIndex *bi, char *key,
                                      BioIndex_IndexFile *index) {
  char * result_string = xcalloc(index->rec_length+1,1);
  Vector * locations     = NULL;
  char * val_ptr       = NULL;
  long   start;
  long   end;
  long   recNum;
  long   i;

  if (BioIndex_binary_search_file(index->handle,
                                  index->rec_length,
                                  key,
                                  index->num_record,
                                  index->header_len,
                                  &recNum,
                                  result_string)) {
    BioIndex_get_match_bounds(bi, key, index, recNum, &start, &end);
    locations = Vector_new();
    for (i=start;i<=end;i++) {
      fseek(index->handle, i*index->rec_length+index->header_len, SEEK_SET);
      fread(result_string, index->rec_length, 1, index->handle);
      result_string[index->rec_length] = '\0';

      if ((val_ptr = BioIndex_get_secondary_value(result_string))!=NULL) {
        Vector_addElement(locations,
                   BioIndex_get_by_primary_key(bi, val_ptr));
      }
    }
  }
                                
  free(result_string);
  return locations;
}


static BioIndex_FileID *BioIndex_FileID_create(char *id, char *name, 
                                               off_t length) {
  BioIndex_FileID *bfid = xcalloc(sizeof(BioIndex_FileID), 1);
  bfid->id = atoi(id);
  bfid->name = StrUtil_copyString(&bfid->name,name,0);
  bfid->length = length;
  return bfid;
}

static void BioIndex_FileID_destroy(BioIndex_FileID *bfid) {
  if (bfid->fp) {
    fclose(bfid->fp);
  }
  free(bfid->name);
  free(bfid);
}
         
static void BioIndex_parse_fileid(BioIndex *bi, char *line) {
  BioIndex_FileID *fileid;
  char **tokens;
  int ntoken;

  StrUtil_tokenizeByDelim(&tokens,&ntoken,line,"\t");

  if (ntoken != 3) {
    Error_write(EOBDA,"BioIndex_parse_fileid",ERR_SEVERE,
                "Format of fileids.dat is incorrect for line %s (ntoken = %d)\n",line,ntoken); 
  }
  fileid = BioIndex_FileID_create(tokens[0],tokens[1],atol(tokens[2]));
  printf("Fileid line = %s %s %s\n",tokens[0],tokens[1],tokens[2]);
  if (bi->fileids == NULL) {
    bi->fileids = Vector_new();
  }
  Vector_addElement(bi->fileids,fileid);
}

static void BioIndex_parse_secondary_namespaces(BioIndex *bi, char *line) {
  int i = 0;
  int nToken = 0;
  char **tokens;

  StrUtil_tokenizeByDelim(&tokens,&nToken,line,"\t");

  bi->num_secondary = nToken;
  bi->secondary_namespaces = xcalloc(sizeof(char *),bi->num_secondary);
  for (i=0;i<nToken;i++) {
    bi->secondary_namespaces[i] = StrUtil_copyString(&(bi->secondary_namespaces[i]),tokens[i],0);
    StrUtil_rmEndSpaces(bi->secondary_namespaces[i]);
    i++;
  }
  bi->num_secondary = nToken;
}

static void BioIndex_parse_primary_namespace(BioIndex *bi, char *line) {
  bi->primary_namespace = StrUtil_copyString(&(bi->primary_namespace),line,0);
  StrUtil_rmEndSpaces(bi->primary_namespace);
}


static void BioIndex_parse_config(BioIndex *bi, char *path) {
  char *config_fname = xcalloc(MAXPATHLEN+1,1);
  FILE *config_fp;
  char line[1024];

  sprintf(config_fname, "%s/config.dat",path);

  printf("File = %s\n",config_fname);
  config_fp = fopen(config_fname, "r");
   
  if (fgets(line,1024,config_fp) == NULL) {
    Error_write(EFGETSNULL,"BioIndex_parse_primary_namespace",ERR_SEVERE,"Failed reading first line from file %s",config_fname);
  }

  bi->fileids = NULL;
  bi->num_secondary = 0;
  while (fgets(line,1024,config_fp)) {
    if (!strncmp(line,"primary_namespace",sizeof("primary_namespace")-1)) {
      BioIndex_parse_primary_namespace(bi,&line[sizeof("primary_namespace")]);
    } else if (!strncmp(line,"fileid_",sizeof("fileid_")-1)) {
      BioIndex_parse_fileid(bi,&line[sizeof("fileid_")]);
    } else if (!strncmp(line,"secondary_namespaces",
                        sizeof("secondary_namespaces")-1)) {
      BioIndex_parse_secondary_namespaces(bi,&line[sizeof("secondary_namespaces")]);
    }
  }
  free(config_fname);
  fclose(config_fp);
}

static void BioIndex_IndexFile_write_header(BioIndex_IndexFile *bif) {
  int i;
  fprintf(bif->handle,"%4d",bif->rec_length);
}

void BioIndex_write_primary(BioIndex *bi, BioIndex_Location *locations, 
                            char *path) {

  char line[1024];
  char format[1024];
  char *primary_fname = xcalloc(MAXPATHLEN+1,1);
  int  i;

  sprintf(primary_fname,"%s/key_%s.key",path,bi->primary_namespace);

  if ((bi->primary_index->handle = fopen(primary_fname,"w")) == NULL) {
    Error_write(EOBDA,"BioIndex_write_primary",ERR_SEVERE,"Failed opening primary index %s\n", primary_fname);
  }

  BioIndex_IndexFile_write_header(bi->primary_index);

  for (i=0; i<bi->primary_index->num_record; i++) {
    BioIndex_Location *loc = (BioIndex_Location *)&locations[i];
    sprintf(format,"%%s\t%%d\t%%ld\t%%d");
    sprintf(line,format,loc->entry, loc->fileIndex,
            loc->start, loc->length);
    sprintf(format,"%%-%d.%ds",bi->primary_index->rec_length,
            bi->primary_index->rec_length);
    fprintf(bi->primary_index->handle, format, line);
  }
  fflush(bi->primary_index->handle);
  free(primary_fname); 
}

typedef struct {
  char *key;
  char *value;
} BioIndex_Secondary_Pair;

typedef struct {
  BioIndex_Index_Definition *definition;
  Vector *ids;
  Vector *final_data;
} BioIndex_Index_Data;

static int BioIndex_secondary_key_compare(const void *a, const void *b) {
  BioIndex_Secondary_Pair *one = *(BioIndex_Secondary_Pair **)a;
  BioIndex_Secondary_Pair *two = *(BioIndex_Secondary_Pair **)b;
  return strcmp((one)->key,(two)->key);
}

static void BioIndex_write_secondaries(BioIndex *bi, Vector *secondary_data,
                                       char *path, int maxValLen) {

  char line[1024];
  char format[1024];
  char *secondary_fname = xcalloc(MAXPATHLEN+1,1);
  int  i;
  int  j;
  int  maxKeyLen;
  int  rec_len;

  if (!Vector_getNumElement(secondary_data)) {
    return;
  }

  bi->secondary_indices = xcalloc(sizeof(BioIndex_IndexFile *),Vector_getNumElement(secondary_data));

  for (i=0;i<Vector_getNumElement(secondary_data);i++) {
    BioIndex_Index_Data *bid = Vector_getElementAt(secondary_data,i);
    sprintf(secondary_fname,"%s/id_%s.index",path,bid->definition->type);

    Vector_sort(bid->final_data, BioIndex_secondary_key_compare);

/*
        for (j=1;j<bid->final_data->len;j++) {
            BioIndex_Secondary_Pair *pair = (BioIndex_Secondary_Pair *)bid->final_data->pdata[j];
            BioIndex_Secondary_Pair *prevpair = (BioIndex_Secondary_Pair *)bid->final_data->pdata[j-1];
            if (!strcmp(pair->key,prevpair->key) && 
                 strcmp(pair->value,prevpair->value)) {
                printf("Multiple values for key %s in secondary\n",pair->key);
                }
            }
*/
    bi->secondary_indices[i] = xcalloc(sizeof(BioIndex_IndexFile),1);
    if ((bi->secondary_indices[i]->handle = fopen(secondary_fname,"w")) == NULL) {
      Error_write(EOBDA,"BioIndex_write_secondaries",ERR_SEVERE,"Failed opening secondary index %s\n", secondary_fname);
    }
    maxKeyLen = 0;
    for (j=0; j<Vector_getNumElement(bid->final_data); j++) {
      BioIndex_Secondary_Pair *pair = (BioIndex_Secondary_Pair *)Vector_getElementAt(bid->final_data,j);
      if (strlen(pair->key) > maxKeyLen) {
        maxKeyLen = strlen(pair->key);
      }
    }
    rec_len = maxKeyLen+maxValLen+1;
    sprintf(format, "%%-%d.%ds",rec_len,rec_len);
    fprintf(bi->secondary_indices[i]->handle,"%4d",rec_len);
    bi->secondary_indices[i]->rec_length = rec_len;
    for (j=0; j<Vector_getNumElement(bid->final_data); j++) {
      BioIndex_Secondary_Pair *pair = (BioIndex_Secondary_Pair *)Vector_getElementAt(bid->final_data,j);
      sprintf(line,"%s\t%s",pair->key, pair->value);
      fprintf(bi->secondary_indices[i]->handle, format, line);
    }
    fflush(bi->secondary_indices[i]->handle);
    free(secondary_fname); 
  }
}

static void BioIndex_write_fileids(BioIndex *bi, FILE *config_fp) {
  int i;

  for (i=0; i<Vector_getNumElement(bi->fileids); i++) {
    BioIndex_FileID *fid = (BioIndex_FileID *)Vector_getElementAt(bi->fileids,i);
    fprintf(config_fp,"fileid_%d\t%s\t" IDFMTSTR "\n",fid->id, fid->name, (long long)fid->length);
  }
}

BioIndex *BioIndex_open_flat(char *path){
  register BioIndex *bi = xcalloc(sizeof(BioIndex), 1);

  bi->path = StrUtil_copyString(&bi->path,path,0);

  /* read config */
  BioIndex_parse_config(bi, path);
    
  /* open index files */
  BioIndex_open_indices(bi, path);

  return bi;
}


int BioIndex_id_compare(const void *a, const void *b) {
  BioIndex_Location *one = (BioIndex_Location *)a;
  BioIndex_Location *two = (BioIndex_Location *)b;
  return strcmp((one)->entry,(two)->entry);
}

BioIndex_Index_Data *BioIndex_Index_Data_create(BioIndex_Index_Definition *def){
  BioIndex_Index_Data *bid = xcalloc(sizeof(BioIndex_Index_Data),1);

  bid->definition = def;
  bid->ids = Vector_new();
  bid->final_data = Vector_new();
  return bid;
}

static BioIndex_Secondary_Pair *BioIndex_Secondary_Pair_create(char *key, char *value) {
  BioIndex_Secondary_Pair *bsp = xcalloc(sizeof(BioIndex_Secondary_Pair),1);

  bsp->key = key;
  bsp->value = value;
  return bsp;
}

static void BioIndex_add_file(BioIndex *bi, char *fname, 
                              BioIndex_Location **locationsP,
                              int index, BioIndex_Index_Definition *primary_def,
                              Vector *secondary_data) {
  FILE *      fpSeq;
/*NOTE DO NOT CHANGE THESE FROM MAXSTRLEN without looking at parsers.c */
  char        line[MAXSTRLEN];
  char        header[MAXSTRLEN];
  long        filePos;
  long        lastPos;
  long        newPos;
  long        startPos;
  long        fileLen;
  char        fileKey[MAXSTRLEN];
  int        first = 1;
  int        i;
  int        j;
  BioIndex_Index_Definition *sbid;
  BioIndex_Index_Data *bid;
  char       *value;

  printf("Match for %s is %s\n",currentPrefix,fname);
  if ((fpSeq = fopen(fname,"r"))  == NULL) {
    Error_write(EOBDA,"BioIndex_add_file",ERR_SEVERE,"Failed opening %s for read",fname);
  }
   
  fseek(fpSeq,0,SEEK_END);
  fileLen = ftell(fpSeq);
  rewind(fpSeq);

  filePos = 0;
  if (fgets(line,MAXSTRLEN,fpSeq) == NULL) {
    Error_write(EFGETSNULL,"BioIndex_add_file",ERR_SEVERE,"Failed reading first line from file %s",fname);
  }
  first = 1;
  lastPos = 0;
   
  while (!feof(fpSeq)) {
    newPos = ftell(fpSeq);
/* For each of the secondaries we need to parse out and store the secondary */
/* keys and wait till the end of the entry (so we known the primary key) */
/* Then add primary key as values into the secondary hashes */
    if (!strncmp(line,primary_def->record_start,primary_def->start_len) || 
        newPos == fileLen) {
      if (newPos == fileLen) {
        printf("Got to end of file\n");
        filePos = fileLen;
      } 
      if (!first) {
        BioIndex_addLocation(bi, locationsP, header, lastPos, filePos-1, index,
                             primary_def->parser);
        for (i=0;i<Vector_getNumElement(secondary_data);i++) {
          BioIndex_Index_Data *bid = Vector_getElementAt(secondary_data,i);
          for (j=Vector_getNumElement(bid->ids)-1;j>=0; j--) {
            Vector_addElement(bid->final_data,
                              BioIndex_Secondary_Pair_create(Vector_getElementAt(bid->ids,j),
                                                             (*locationsP)[bi->primary_index->num_record-1].entry));
            Vector_removeElementAt(bid->ids,j);
          }
        } 
      } else {
        first = 0;
      }
      lastPos = filePos;
    }
    if (!strncmp(line,primary_def->line_prefix,primary_def->prefix_len)) {
      strcpy(header, line);
    }
    for (i=0;i<Vector_getNumElement(secondary_data);i++) {
      BioIndex_Index_Data *bid = Vector_getElementAt(secondary_data,i);
      sbid = bid->definition;
      if (!strncmp(line,sbid->line_prefix,sbid->prefix_len)) {
        if (sbid->multiParser) {
          char *id;
          char *pos = NULL;
          char *current = line;
          while ((id = sbid->multiParser(current,&pos))) {
            Vector_addElement(bid->ids,id);
            current = pos;
          }
        } else {
          Vector_addElement(bid->ids,sbid->parser(line));
        }
      }
    }
    filePos = newPos;
    // Don't really need to check for return value here as going to do an feof check at the top of the loop - 
    // get the return value to satify the compiler
    char *ret = fgets(line,MAXSTRLEN,fpSeq);
  }
}         


int selectPrefix(const struct dirent *DirEnt) {
  int lenPref = strlen(currentPrefix);

  if (!strncmp(DirEnt->d_name,currentPrefix,lenPref)) {
    if (currentPrefix[lenPref-1] == '-' && 
        areintquiet(DirEnt->d_name+lenPref)) {
      return 1;
    } else if (strlen(DirEnt->d_name) == lenPref) {
      return 1;
    }
  }
  return 0;
}   

static void BioIndex_write_config(BioIndex *bi, char *path, char *type,
                                  BioIndex_Index_Definition *primary_def,
                                  Vector *secondary_defs) {
  FILE *fp;
  char fileName[MAXPATHLEN];
  int i;

  sprintf(fileName,"%s/%s",path,"config.dat");
  if ((fp = fopen(fileName,"wb")) == NULL) {
    Error_write(EOBDA,"BioIndex_write_config",ERR_SEVERE,"Failed opening %s for write\n",fileName);
  }

  fprintf(fp,"index\t%s\n",type);
  fprintf(fp,"primary_namespace\t%s\n",primary_def->type);
  if (bi->num_secondary) {
    fprintf(fp,"secondary_namespaces");
    for (i=0;i<bi->num_secondary;i++) {
      BioIndex_Index_Definition *bid = Vector_getElementAt(secondary_defs,i);
      fprintf(fp,"\t%s",bid->type);
    }
    fprintf(fp,"\n");
  }
  BioIndex_write_fileids(bi, fp);
  fclose(fp);
}

BioIndex *BioIndex_generate_flat(char *path, char *seq_path, 
                                 char *select, char *format,
                                 BioIndex_Index_Definition *primary_def,
                                 Vector *secondary_defs) {
  char  fileName[1024];
  FILE  *fp;
  struct dirent **FileNames; 
  int    NFile;
  char   FullFName[1024];  
  int    i;
  BioIndex_Location *locations = xcalloc(sizeof(BioIndex_Location),1000000);
  register BioIndex *bi = xcalloc(sizeof(BioIndex), 1);
  char numString[1024];
  long maxStart = 0;
  int maxLength = 0;
  int maxIdLen  = 0;
  int maxIndex  = 0;
  struct stat stats;
  BioIndex_FileID *fileid;
  Vector  *secondary_data = Vector_new();
  int fieldpad  = 0;
    
  bi->primary_namespace = StrUtil_copyString(&bi->primary_namespace,primary_def->type,0);

  currentPrefix = select;
  if ((NFile=scandir(seq_path,&FileNames,selectPrefix,alphasort))== -1) {
    Error_write(EOBDA,"BioIndex_generate_flat",ERR_SEVERE,"Couldn't open sequence directory %s",seq_path);
  } 

  if (NFile == 0) {
    Error_write(EOBDA,"BioIndex_generate_flat",ERR_SEVERE,"Didn't find any files matching %s\n",currentPrefix);
  }

  for (i=0;i<Vector_getNumElement(secondary_defs);i++) {
    Vector_addElement(secondary_data, BioIndex_Index_Data_create(Vector_getElementAt(secondary_defs,i)));
  }

  bi->fileids = Vector_new();
  bi->primary_index = xcalloc(sizeof(BioIndex_IndexFile), 1);
  for (i=0; i<NFile; i++) {
    getFullFName(FullFName, seq_path, FileNames[i]->d_name);
    if (stat(FullFName,&stats) == 0) {
      sprintf(numString,"%d",i);
      printf("Size of %s = " IDFMTSTR "\n",FullFName,(long long)stats.st_size);
      fileid = BioIndex_FileID_create(numString,FullFName,
                                      stats.st_size);
      Vector_addElement(bi->fileids,fileid);
      BioIndex_add_file(bi, FullFName, &locations, i, primary_def,
                        secondary_data);
    } else {
      perror(FullFName);
    }
  }   

/* sort */
  qsort(locations,bi->primary_index->num_record, sizeof(BioIndex_Location),
        BioIndex_id_compare);
 
  for (i=1; i<bi->primary_index->num_record; i++) {
    BioIndex_Location *loc = (BioIndex_Location *)&locations[i];
    BioIndex_Location *prev_loc = (BioIndex_Location *)&locations[i-1];

    if (!strcmp(loc->entry,prev_loc->entry)) {
      Error_write(EOBDA,"BioIndex_generate_flat",ERR_SEVERE,"Primary key not unique for %s (%ld) %s (%ld)\n",
              loc->entry,loc->start,prev_loc->entry,prev_loc->start);
    }          
  }
/* Find maximal widths for fields */
  
  bi->primary_index->num_field = 4;
  bi->primary_index->field_properties = xcalloc(sizeof(BioIndex_FieldProperties),4);

  for (i=0; i<bi->primary_index->num_record; i++) {
    BioIndex_Location *loc = (BioIndex_Location *)&locations[i];
    if (strlen(loc->entry) > maxIdLen) {
      maxIdLen = strlen(loc->entry);
    }
    if (loc->start > maxStart) {
      maxStart = loc->start;
    }
    if (loc->length > maxLength) {
      maxLength = loc->length;
    }
    if (loc->fileIndex > maxIndex) {
      maxIndex = loc->fileIndex;
    }
  }


  bi->primary_index->rec_length = 
       (bi->primary_index->field_properties[0].width = maxIdLen + fieldpad);
  bi->primary_index->rec_length += 
       (bi->primary_index->field_properties[1].width = sprintf(numString,"%d",maxIndex)+fieldpad);
  bi->primary_index->rec_length += 
       (bi->primary_index->field_properties[2].width = sprintf(numString,"%ld",maxStart)+fieldpad);
  bi->primary_index->rec_length += 
       (bi->primary_index->field_properties[3].width = sprintf(numString,"%d",maxLength)+fieldpad);
  bi->primary_index->rec_length += 3; /* for tabs */

  bi->primary_index->field_properties[0].type = 's';
  bi->primary_index->field_properties[1].type = 'd';
  bi->primary_index->field_properties[2].type = 'l';
  bi->primary_index->field_properties[3].type = 'd';

/* Write out records */
  BioIndex_write_primary(bi, locations, path);
  BioIndex_write_secondaries(bi, secondary_data, path, maxIdLen);
  bi->num_secondary = Vector_getNumElement(secondary_data);

  BioIndex_write_config(bi, path, "flat/1", primary_def, secondary_defs);

  free(locations);

  return bi;
}


void BioIndex_close(BioIndex *bi){
  int i;

  free(bi->path);

  if (bi->primary_index->handle) {
    fclose(bi->primary_index->handle);
  }

  for (i=0;i<Vector_getNumElement(bi->fileids); i++) {
    BioIndex_FileID_destroy(Vector_getElementAt(bi->fileids,i));
  }
  Vector_free(bi->fileids);

  free(bi);
  return;
}

int areintquiet(const char *str) {
  int length = 0;
  int i;
 
  if ((length = strlen(str))==0) {
    return 0;
  }
  for (i=0;i<length;i++) {
    if (!isdigit(str[i])&&str[i]!='-') {
      return 0;
    }
  }
  return 1;
}


BioIndex_Index_Definition *BioIndex_Index_Definition_create(char *record_start,
                                              char *line_prefix,
                                              BioIndex_Parser_Func parser,
                                              BioIndex_MultiParser_Func multi_parser,
                                              char *type){
  register BioIndex_Index_Definition *def = xcalloc(sizeof(BioIndex_Index_Definition),1);

  def->record_start = StrUtil_copyString(&def->record_start,record_start,0);
  def->line_prefix  = StrUtil_copyString(&def->line_prefix,line_prefix,0);
  def->prefix_len   = strlen(line_prefix);
  def->start_len    = strlen(record_start);
  def->parser       = parser;
  def->multiParser  = multi_parser;
  def->type         = StrUtil_copyString(&def->type,type,0);

  return def;
}
