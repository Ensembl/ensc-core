#include "AssemblyMapperAdaptor.h"
#include "BaseAdaptor.h"
#include "MysqlUtil.h"
#include "AssemblyMapper.h"
#include "GenomicRange.h"

#include "StatementHandle.h"
#include "ResultRow.h"



AssemblyMapperAdaptor *AssemblyMapperAdaptor_new(DBAdaptor *dba) {
  AssemblyMapperAdaptor *ama;

  if ((ama = (AssemblyMapperAdaptor *)calloc(1,sizeof(AssemblyMapperAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for AssemblyMapperAdaptor\n");
    return NULL;
  }
  BaseAdaptor_init((BaseAdaptor *)ama, dba, ASSEMBLYMAPPER_ADAPTOR);
  ama->typeCache = StringHash_new(STRINGHASH_SMALL);

  return ama;
}

AssemblyMapper *AssemblyMapperAdaptor_fetchByType(AssemblyMapperAdaptor *ama, char *type) {

  if( !StringHash_contains(ama->typeCache,type)) {
    StringHash_add(ama->typeCache,type,AssemblyMapper_new(ama,type));
  }

  return StringHash_getValue(ama->typeCache,type);
}


void AssemblyMapperAdaptor_registerRegion(AssemblyMapperAdaptor *ama, 
                                          AssemblyMapper *assMapper,
                                          char *assemblyType,
                                          IDType chrId,
                                          int chrStart,
                                          int chrEnd) {
  char qStr[512];
  StatementHandle *sth;
  ResultRow *row;

  sprintf(qStr,
    "select"
    "    ass.contig_start,"
    "    ass.contig_end,"
    "    ass.contig_id,"
    "    ass.contig_ori,"
    "    chr.chromosome_id,"
    "    ass.chr_start,"
    "    ass.chr_end"
    " from"
    "    assembly ass,"
    "    chromosome chr"
    " where"
    "    chr.chromosome_id = "
    IDFMTSTR
    " and "
    "    ass.chromosome_id = chr.chromosome_id and"
    "    %d <= ass.chr_end  and"
    "    %d >= ass.chr_start  and"
    "    ass.type = '%s'",
    chrId,chrStart,chrEnd,assemblyType);

  sth = ama->prepare((BaseAdaptor *)ama,qStr,strlen(qStr));
  sth->execute(sth);

  while (row = sth->fetchRow(sth)) {
    IDType contigId = row->getLongLongAt(row,2);

    if (!AssemblyMapper_haveRegisteredContig(assMapper,contigId)) {
      Mapper *mapper = AssemblyMapper_getMapper(assMapper);

      AssemblyMapper_registerContig(assMapper,contigId);
      
      Mapper_addMapCoordinates(mapper,contigId,
                               row->getIntAt(row,0),
                               row->getIntAt(row,1),
                               row->getIntAt(row,3),
                               row->getLongLongAt(row,4),
                               row->getIntAt(row,5),
                               row->getIntAt(row,6));
    }
  }
  sth->finish(sth);

  return;
}


GenomicRange *AssemblyMapperAdaptor_registerContig(AssemblyMapperAdaptor *ama, 
                                          AssemblyMapper *assMapper,
                                          char *assemblyType,
                                          IDType contigId) {
  char qStr[512];
  StatementHandle *sth;
  ResultRow *row;
  int extraRows = 0;
  GenomicRange *range;

  sprintf(qStr,
    "select"
    "  c.name,"
    "  a.chr_start,"
    "  a.chr_end,"
    "  c.chromosome_id"
    " from"
    "  assembly a,"
    "  chromosome c"
    " where"
    "   contig_id = "
    IDFMTSTR
    " and"
    "   type = '%s' and"
    "   c.chromosome_id = a.chromosome_id",
    contigId, assemblyType);

  sth = ama->prepare((BaseAdaptor *)ama,qStr,strlen(qStr));
  sth->execute(sth);

  row = sth->fetchRow(sth);

  if (row == NULL) {
    sth->finish(sth);
    return NULL;
  }

  range = GenomicRange_new();

  GenomicRange_setChrName(range,row->getStringAt(row,0));
  GenomicRange_setChrStart(range,row->getIntAt(row,1));
  GenomicRange_setChrEnd(range,row->getIntAt(row,2));
  GenomicRange_setChrId(range,row->getLongAt(row,3));

  while (row = sth->fetchRow(sth)) {
    extraRows++;
  }

  sth->finish(sth);

  if (extraRows) {
    fprintf(stderr,"WARNING: Contig %d is ambiguous in assembly type %s\n",contigId, assemblyType);
    GenomicRange_free(range);
    return NULL;
  }

  return range;
}
