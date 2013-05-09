#include "SupportingFeatureAdaptor.h"

#include "DNAAlignFeatureAdaptor.h"
#include "ProteinAlignFeatureAdaptor.h"
#include "BaseAlignFeature.h"
#include "DBAdaptor.h"
#include "BaseFeatureAdaptor.h"

#include <string.h>

SupportingFeatureAdaptor *SupportingFeatureAdaptor_new(DBAdaptor *dba) {
  SupportingFeatureAdaptor *sfa;

  if ((sfa = (SupportingFeatureAdaptor *)calloc(1,sizeof(SupportingFeatureAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for SupportingFeatureAdaptor\n");
    return NULL;
  }
  BaseAdaptor_init((BaseAdaptor *)sfa, dba, SUPPORTINGFEATURE_ADAPTOR);

  return sfa;
}



Vector *SupportingFeatureAdaptor_fetchAllByExon(SupportingFeatureAdaptor *sfa, Exon *exon) {
  StatementHandle *sth;
  char qStr[512];
  Vector *out = Vector_new();
  ResultRow *row;
  DNAAlignFeatureAdaptor *dafa;
  ProteinAlignFeatureAdaptor *pafa;

  // throw if this is a sticky exon
  if (exon->objectType == CLASS_STICKYEXON) {
    fprintf(stderr,"Expected Exon but got StickyExon. Call get_all_component_Exons first\n");
    exit(1);
  }

  if (!Exon_getDbID(exon)) {
    fprintf(stderr,"WARNING: exon has no dbID can't fetch evidence from db " 
                   "no relationship exists\n");
    return emptyVector;
  }

  sprintf(qStr,"SELECT sf.feature_type, sf.feature_id "
               "FROM   supporting_feature sf "
               "WHERE  exon_id = " IDFMTSTR, Exon_getDbID(exon));

  sth = sfa->prepare((BaseAdaptor *)sfa, qStr, strlen(qStr));

  sth->execute(sth);
  
  pafa = DBAdaptor_getProteinAlignFeatureAdaptor(sfa->dba);
  dafa = DBAdaptor_getDNAAlignFeatureAdaptor(sfa->dba);

  
  while ((row = sth->fetchRow(sth))) {      
    BaseAlignFeature *baf = NULL;
    char *type = row->getStringAt(row,0);

// baf is HACK HACK HACK
    if (!strcmp(type,"protein_align_feature")) {
      baf = (BaseAlignFeature *)ProteinAlignFeatureAdaptor_fetchByDbID(pafa, row->getLongLongAt(row,1));
    } else if (!strcmp(type,"dna_align_feature")) {
      baf = (BaseAlignFeature *)DNAAlignFeatureAdaptor_fetchByDbID(dafa, row->getLongLongAt(row,1));
    } else {
      fprintf(stderr,"Error: Unknown feature type [%s]\n",type);
      exit(1);
    }

    if (Exon_getContig(exon)->objectType == CLASS_SLICE) {
      //tranform to slice coords
      BaseAlignFeature_transformToSlice(baf, (Slice *)Exon_getContig(exon));
    } else {
      //we might have to convert the features coordinate system
      if (BaseContig_getDbID(BaseAlignFeature_getContig(baf)) != BaseContig_getDbID(Exon_getContig(exon))) {
        // fprintf(stderr,"Error: Free for different contig skip not implemented\n");
        BaseAlignFeature_free(baf);
        continue;
      }
    }

    Vector_addElement(out,baf);
  }

  sth->finish(sth);
  return out;
}
