#include "ChromosomeAdaptor.h"
#include "BaseAdaptor.h"
#include "MysqlUtil.h"

#include "StatementHandle.h"
#include "ResultRow.h"


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

Chromosome *ChromosomeAdaptor_fetchByDbID(ChromosomeAdaptor *ca, IDType dbID) {
  Chromosome *chromosome;
  char qStr[256];
  StatementHandle *sth;
  ResultRow *row;

  if (IDHash_contains(ca->chrCache,dbID)) {

    chromosome = IDHash_getValue(ca->chrCache, dbID);

  } else {
    sprintf(qStr,"SELECT chromosome_id, name, length"
      " FROM chromosome"
      " WHERE  chromosome_id = "
      IDFMTSTR, dbID);
  
    sth = ca->prepare((BaseAdaptor *)ca,qStr,strlen(qStr));
    sth->execute(sth);
  
    row = sth->fetchRow(sth);
    if( row == NULL ) {
      sth->finish(sth);
      return NULL;
    }
  
    chromosome = ChromosomeAdaptor_chromosomeFromRow(ca, row);
    sth->finish(sth);
  }

  return chromosome;
}

Chromosome *ChromosomeAdaptor_fetchByChrName(ChromosomeAdaptor *ca, char *chrName) {
  Chromosome *chromosome;
  char qStr[256];
  StatementHandle *sth;
  ResultRow *row;

  if (StringHash_contains(ca->chrNameCache,chrName)) {

    chromosome = StringHash_getValue(ca->chrNameCache, chrName);

  } else {
    sprintf(qStr,"SELECT chromosome_id, name, length"
      " FROM chromosome"
      " WHERE  name = '%s'", chrName);
  
    sth = ca->prepare((BaseAdaptor *)ca,qStr,strlen(qStr));
    sth->execute(sth);
  
    row = sth->fetchRow(sth);
    if( row == NULL ) {
      sth->finish(sth);
      fprintf(stderr, "ERROR: Do not recognise chromosome %s\n",chrName);
      exit(1);
    }
  
    chromosome = ChromosomeAdaptor_chromosomeFromRow(ca, row);
    sth->finish(sth);
  }

  return chromosome;
}

Chromosome *ChromosomeAdaptor_chromosomeFromRow(ChromosomeAdaptor *ca, 
                                                ResultRow *row) {
  Chromosome *chromosome = Chromosome_new();

  Chromosome_setAdaptor(chromosome,(BaseAdaptor *)ca);
  Chromosome_setDbID(chromosome,row->getLongLongAt(row,0));
  Chromosome_setName(chromosome,row->getStringAt(row,1));
  Chromosome_setLength(chromosome,row->getIntAt(row,2));

  IDHash_add(ca->chrCache, Chromosome_getDbID(chromosome), chromosome);
  StringHash_add(ca->chrNameCache, Chromosome_getName(chromosome), chromosome);
  return chromosome;
}
