#include "ChromosomeAdaptor.h"
#include "BaseAdaptor.h"
#include "MysqlUtil.h"


ChromosomeAdaptor *ChromosomeAdaptor_new(DBAdaptor *dba) {
  ChromosomeAdaptor *ca;

  if ((ca = (ChromosomeAdaptor *)calloc(1,sizeof(ChromosomeAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for ChromosomeAdaptor\n");
    return NULL;
  }
  BaseAdaptor_init((BaseAdaptor *)ca, dba, CHROMOSOME_ADAPTOR);

  ca->chrCache = IDHash_new(IDHASH_SMALL);
  ca->chrNameCache = StringHash_new(STRINGHASH_SMALL);

  return ca;
}

Chromosome *ChromosomeAdaptor_fetchByDbID(ChromosomeAdaptor *ca, long dbID) {
  Chromosome *chromosome;
  char qStr[256];
  MYSQL_RES *results;
  MYSQL_ROW row;

  if (IDHash_contains(ca->chrCache,dbID)) {

    chromosome = IDHash_getValue(ca->chrCache, dbID);

  } else {
    sprintf(qStr,"SELECT chromosome_id, name, length"
      " FROM chromosome"
      " WHERE  chromosome_id = %d", dbID);
  
    results = ca->prepare((BaseAdaptor *)ca,qStr,strlen(qStr));
  
    row = mysql_fetch_row(results);
    if( row == NULL ) {
      return NULL;
    }
  
    chromosome = ChromosomeAdaptor_chromosomeFromRow(ca, row);
  }

  return chromosome;
}

Chromosome *ChromosomeAdaptor_fetchByChrName(ChromosomeAdaptor *ca, char *chrName) {
  Chromosome *chromosome;
  char qStr[256];
  MYSQL_RES *results;
  MYSQL_ROW row;

  if (StringHash_contains(ca->chrNameCache,chrName)) {

    printf("Using chromosome name cache\n");
    chromosome = StringHash_getValue(ca->chrNameCache, chrName);

  } else {
    sprintf(qStr,"SELECT chromosome_id, name, length"
      " FROM chromosome"
      " WHERE  name = '%s'", chrName);
  
    results = ca->prepare((BaseAdaptor *)ca,qStr,strlen(qStr));
  
    row = mysql_fetch_row(results);
    if( row == NULL ) {
      fprintf(stderr, "ERROR: Do not recognise chromosome %s\n",chrName);
      exit(1);
    }
  
    chromosome = ChromosomeAdaptor_chromosomeFromRow(ca, row);
  }

  return chromosome;
}

Chromosome *ChromosomeAdaptor_chromosomeFromRow(ChromosomeAdaptor *ca, 
                                                MYSQL_ROW row) {
  Chromosome *chromosome = Chromosome_new();

  Chromosome_setAdaptor(chromosome,(BaseAdaptor *)ca);
  Chromosome_setDbID(chromosome,MysqlUtil_getLong(row,0));
  Chromosome_setName(chromosome,MysqlUtil_getString(row,1));
  Chromosome_setLength(chromosome,MysqlUtil_getInt(row,2));

  IDHash_add(ca->chrCache, Chromosome_getDbID(chromosome), chromosome);
  StringHash_add(ca->chrNameCache, Chromosome_getName(chromosome), chromosome);
  return chromosome;
}
