#include "SyntenyAdaptor.h"

#include "StrUtil.h"

SyntenyAdaptor *SyntenyAdaptor_new(ComparaDBAdaptor *dba) {
  SyntenyAdaptor *sa;

  if ((sa = (SyntenyAdaptor *)calloc(1,sizeof(SyntenyAdaptor))) == NULL) {
    fprintf(stderr,"Error: Failed allocating sa\n");
    exit(1);
  }
  BaseComparaAdaptor_init((BaseComparaAdaptor *)sa, dba, SYNTENY_ADAPTOR);

  return sa;
}

void SyntenyAdaptor_setSpecies(SyntenyAdaptor *sa, char *species1, char *species2) {

    sa->speciesMain      = StrUtil_copyString(&sa->speciesMain,species1,0);
    StrUtil_strReplChr(sa->speciesMain,'_',' ');
    sa->speciesSecondary = StrUtil_copyString(&sa->speciesSecondary,species2,0);
    StrUtil_strReplChr(sa->speciesSecondary,'_',' ');
}

// if chr = NULL return all synteny pairs
Vector *SyntenyAdaptor_getSyntenyForChromosome(SyntenyAdaptor *sa, char *chr, int *startP, int *endP) {
  Vector *data;
  char extraSql[128] = "";
  char qStr[2048];
  ResultRow *row;
  StatementHandle *sth;

  if (!strcmp(SyntenyAdaptor_getSpeciesMain(sa), SyntenyAdaptor_getSpeciesSecondary(sa))) { 
    return emptyVector;
  }


  if (chr) {
    sprintf(extraSql, " and df.name = '%s'",chr);
    if (startP) {
      char tmpStr[128];
      sprintf(tmpStr," and dfr.seq_start <= %d and dfr.seq_end >= %d",*endP,*startP);
      strcat(extraSql,tmpStr);
    }
  }

  sprintf(qStr,
        "select sr.synteny_region_id,"
        "       df.dnafrag_type as core_type,  df.name as core_name,"
        "       dfr.seq_start as core_start,   dfr.seq_end as core_end,"
        "       df_h.dnafrag_type as hit_type, df_h.name as hit_name,"
        "       dfr_h.seq_start as hit_start,  dfr_h.seq_end as hit_end,"
        "       sr.rel_orientation"
        "  from dnafrag as df,         dnafrag as df_h,"
        "       dnafrag_region as dfr, dnafrag_region as dfr_h,"
        "       genome_db as gd,       genome_db as gd_h,"
        "       synteny_region as sr"
        " where gd.name = '%s'   and gd.genome_db_id   = df.genome_db_id and"
        "       gd_h.name = '%s' and gd_h.genome_db_id = df_h.genome_db_id and"
        "       df.dnafrag_id   = dfr.dnafrag_id and"
        "       df_h.dnafrag_id = dfr_h.dnafrag_id and"
        "       dfr.synteny_region_id   = sr.synteny_region_id and"
        "       dfr_h.synteny_region_id = sr.synteny_region_id %s"
        " order by df.name, dfr.seq_start", 
          SyntenyAdaptor_getSpeciesMain(sa),
          SyntenyAdaptor_getSpeciesSecondary(sa),
          extraSql);

  sth = sa->prepare((BaseAdaptor *)sa, qStr, strlen(qStr));
  sth->execute(sth);

  data = Vector_new();

  while ((row = sth->fetchRow(sth))) {
    SyntenyRegion *sr = SyntenyRegion_new();

    SyntenyRegion_setDbID(sr, row->getLongLongAt(row,0));
    SyntenyRegion_setSeqType(sr, row->getStringAt(row,1));
    SyntenyRegion_setChrName(sr, row->getStringAt(row,2));
    SyntenyRegion_setChrStart(sr, row->getIntAt(row,3));
    SyntenyRegion_setChrEnd(sr, row->getIntAt(row,4));
    if (startP && endP) {
      SyntenyRegion_setStart(sr, row->getIntAt(row,3)- *startP +1);
      SyntenyRegion_setEnd(sr, row->getIntAt(row,4)- *startP +1);
    } else {
      SyntenyRegion_setStart(sr, row->getIntAt(row,3));
      SyntenyRegion_setEnd(sr, row->getIntAt(row,4));
    }
    SyntenyRegion_setHitSeqType(sr, row->getStringAt(row,5));
    SyntenyRegion_setHitChrName(sr, row->getStringAt(row,6));
    SyntenyRegion_setHitChrStart(sr, row->getIntAt(row,7));
    SyntenyRegion_setHitChrEnd(sr, row->getIntAt(row,8));
    SyntenyRegion_setRelOri(sr, row->getIntAt(row,9));

    Vector_addElement(data, sr);
  }

  sth->finish(sth);
  return data;
}
