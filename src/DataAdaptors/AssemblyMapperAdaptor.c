#include "AssemblyMapperAdaptor.h"
#include "BaseAdaptor.h"
#include "MysqlUtil.h"
#include "AssemblyMapper.h"
#include "GenomicRange.h"



AssemblyMapperAdaptor *AssemblyMapperAdaptor_new(DBAdaptor *dba) {
  AssemblyMapperAdaptor *ama;

  if ((ama = (AssemblyMapperAdaptor *)calloc(1,sizeof(AssemblyMapperAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for AssemblyMapperAdaptor\n");
    return NULL;
  }
  BaseAdaptor_init((BaseAdaptor *)ama, dba, ASSEMBLYMAPPER_ADAPTOR);

  return ama;
}

void AssemblyMapperAdaptor_registerRegion(AssemblyMapperAdaptor *ama, 
                                          AssemblyMapper *assMapper,
                                          char *assemblyType,
                                          long chrId,
                                          int chrStart,
                                          int chrEnd) {
  char qStr[512];
  MYSQL_RES *results;
  MYSQL_ROW row;

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
    "    chr.chromosome_id= '%d' and"
    "    ass.chromosome_id = chr.chromosome_id and"
    "    %d <= ass.chr_end  and"
    "    %d >= ass.chr_start  and"
    "    ass.type = '%s'",
    chrId,chrStart,chrEnd,assemblyType);

  results = ama->prepare((BaseAdaptor *)ama,qStr,strlen(qStr));

  while (row = mysql_fetch_row(results)) {
    long contigId = MysqlUtil_getLong(row,2);

    if (!AssemblyMapper_haveRegisteredContig(assMapper,contigId)) {
      Mapper *mapper = AssemblyMapper_getMapper(assMapper);

      AssemblyMapper_registerContig(assMapper,contigId);
      
      Mapper_addMapCoordinates(mapper,contigId,
                               MysqlUtil_getInt(row,0),
                               MysqlUtil_getInt(row,1),
                               MysqlUtil_getInt(row,3),
                               MysqlUtil_getLong(row,3),
                               MysqlUtil_getInt(row,4),
                               MysqlUtil_getInt(row,5));
    }
  }

  return;
}


GenomicRange *AssemblyMapperAdaptor_registerContig(AssemblyMapperAdaptor *ama, 
                                          AssemblyMapper *assMapper,
                                          char *assemblyType,
                                          long contigId) {
  char qStr[512];
  MYSQL_RES *results;
  MYSQL_ROW row;
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
    "   contig_id = %d and"
    "   type = '%s' and"
    "   c.chromosome_id = a.chromosome_id",
    contigId, assemblyType);

  results = ama->prepare((BaseAdaptor *)ama,qStr,strlen(qStr));

  row = mysql_fetch_row(results);

  if (row == NULL) {
    return NULL;
  }

  range = GenomicRange_new();

  GenomicRange_setChrName(range,MysqlUtil_getString(row,0));
  GenomicRange_setChrStart(range,MysqlUtil_getInt(row,1));
  GenomicRange_setChrEnd(range,MysqlUtil_getInt(row,2));
  GenomicRange_setChrId(range,MysqlUtil_getLong(row,3));

  while (row = mysql_fetch_row(results)) {
    extraRows++;
  }

  if (extraRows) {
    fprintf(stderr,"WARNING: Contig %d is ambiguous in assembly type %s\n",contigId, assemblyType);
    GenomicRange_free(range);
    return NULL;
  }

  return range;
}
