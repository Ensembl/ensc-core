#include "SyntenyAdaptor.h"


void SyntenyAdaptor_setSpecies(SyntenyAdaptor *sa, char *species1, char *species2) {

    sa->speciesMain      = StrUtil_copyString(&sa->speciesMain,species1,0);
    StrUtil_strReplChr(sa->speciesMain,'_',' ');
    sa->speciesSecondary = StrUtil_copyString(&sa->speciesSecondary,species2,0);
    StrUtil_strReplChr(sa->speciesSecondary,'_',' ');
}

// if chr = NULL return all synteny pairs
Vector *SyntenyApdaptor_getSyntenyForChromosome(SyntenyAdaptor *sa, char *chromosome, int *startP, int *endP) {
  Vector *data;
  char extraSql[128] = "";
  char qStr[2048];

  if (!strcmp(SyntenyAdaptor_getSpeciesMain(sa), SyntenyAdaptor_getSpeciesSecondary(sa))) { 
    return emptyVector;
  }


  if (chr) {
    sprintf(extraSql, " and df.name = '%s'",chr);
    if (*startP) {
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
          DNAAlignFeatureAdaptor_getSpeciesMain(dafa),
          DNAAlignFeatureAdaptor_getSpeciesSecondary(dafa),
          extraSql);

  sth = sa->prepare((BaseAdaptor *)sa, qStr, strlen(qStr));
  sth->execute(sth);

  while ((row = sth->fetchRow(sth))) {
    SyntenyRegion *sr = SyntenyRegion_new();

            'synteny_id'    => $Q->[0],
            'seq_type'      => $Q->[1],
            'chr_name'      => $Q->[2],
            'start'         => $Q->[3],
            'end'           => $Q->[4],
            'chr_start'     => $Q->[3],
            'chr_end'       => $Q->[4],
            'start'         => $Q->[3]-$start+1, 
            'end'           => $Q->[4]-$start+1,
            'hit_seq_type'  => $Q->[5],
            'hit_chr_name'  => $Q->[6],
            'hit_chr_start' => $Q->[7],
            'hit_chr_end'   => $Q->[8],
            'rel_ori'       => $Q->[9]
    Vector_addElement(data, sr);
  }
  return data;
}
