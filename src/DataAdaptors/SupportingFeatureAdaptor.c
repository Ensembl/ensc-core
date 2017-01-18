/*
 * Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
 * Copyright [2016-2017] EMBL-European Bioinformatics Institute
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *      http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/*
=head1 NAME

Bio::EnsEMBL::DBSQL::SupportingFeatureAdaptor - Retrieves supporting
features from the database.
*/


#include "SupportingFeatureAdaptor.h"

#include "DNAAlignFeatureAdaptor.h"
#include "ProteinAlignFeatureAdaptor.h"
#include "BaseAlignFeature.h"
#include "DNAAlignFeature.h"
#include "DNAPepAlignFeature.h"
#include "DBAdaptor.h"
#include "BaseFeatureAdaptor.h"
#include "SliceAdaptor.h"
#include "SeqFeature.h"

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


/*
=head2 fetch_all_by_Exon

  Arg [1]    : Bio::EnsEMBL::Exon $exon 
               The exon to fetch supporting features.
  Example    : @sfs =
                @{ $supporting_feat_adaptor->fetch_all_by_Exon($exon) };
  Description: Retrieves supporting features (evidence) for a given
               exon.
  Returntype : List of Bio::EnsEMBL::BaseAlignFeatures in the same
               coordinate system as the $exon argument
  Exceptions : Warning if $exon is not in the database (i.e. dbID
               not defined).
               Throw if a retrieved supporting feature is of unknown
               type.
  Caller     : Bio::EnsEMBL::Exon
  Status     : Stable

=cut
*/
Vector *SupportingFeatureAdaptor_fetchAllByExon(SupportingFeatureAdaptor *sfa, Exon *exon) {
  StatementHandle *sth;
  char qStr[512];
  ResultRow *row;
  DNAAlignFeatureAdaptor *dafa;
  ProteinAlignFeatureAdaptor *pafa;


  if (!Exon_getDbID(exon)) {
    fprintf(stderr,"WARNING: exon has no dbID can't fetch evidence from db "
                   "no relationship exists\n");
    return Vector_new();
  }

  Vector *out = Vector_new();
  sprintf(qStr,"SELECT sf.feature_type, sf.feature_id "
               "FROM   supporting_feature sf "
               "WHERE  exon_id = " IDFMTSTR, Exon_getDbID(exon));

  sth = sfa->prepare((BaseAdaptor *)sfa, qStr, strlen(qStr));

  sth->execute(sth);

  pafa = DBAdaptor_getProteinAlignFeatureAdaptor(sfa->dba);
  dafa = DBAdaptor_getDNAAlignFeatureAdaptor(sfa->dba);

  while ((row = sth->fetchRow(sth))) {
    SeqFeature *sf = NULL;
    char *type = row->getStringAt(row,0);
    IDType dbId = row->getLongLongAt(row,1);

// sf is HACK HACK HACK
    if (!strcmp(type,"protein_align_feature")) {
      sf = (SeqFeature *)ProteinAlignFeatureAdaptor_fetchByDbID(pafa, dbId);
    } else if (!strcmp(type,"dna_align_feature")) {
      sf = (SeqFeature *)DNAAlignFeatureAdaptor_fetchByDbID(dafa, dbId);
    } else {
      fprintf(stderr,"Error: Unknown feature type [%s]\n",type);
      exit(1);
    }

    if (sf == NULL) {
      fprintf(stderr,"Warning: Supporting feature %s "IDFMTSTR" does not exist in DB\n", type, dbId);
    } else {
      SeqFeature *newSf = SeqFeature_transfer(sf, (Slice *)Exon_getSlice(exon));
      // NIY: Free pretranferred one??
      Vector_addElement(out, newSf);
    }
  }

  sth->finish(sth);
  return out;
}


Vector *SupportingFeatureAdaptor_fetchAllByExonList(SupportingFeatureAdaptor *sfa, Vector *exons, Slice *slice) {
  int i;

  // associate exon identifiers with transcripts
  IDHash *exHash = IDHash_new(IDHASH_MEDIUM);
  for (i=0; i<Vector_getNumElement(exons); i++) {
    Exon *e  = Vector_getElementAt(exons, i);
    if ( ! IDHash_contains(exHash, Exon_getDbID(e))) {
      Vector *exVec = Vector_new();
      IDHash_add(exHash, Exon_getDbID(e), exVec);
    }
    // I thought exons would always have been pruned to one exon for each db id, but it appears that if transcripts are 
    // transferred (can happen in gene adaptor fetch by slice), then duplicate shallow copies are made.
    Vector *exVec = IDHash_getValue(exHash, Exon_getDbID(e));
    Vector_addElement(exVec, e);
/*
    } else {
      Exon *compEx = IDHash_getValue(exHash, Exon_getDbID(e));
      if (compEx != e) {
        fprintf(stderr,"ERROR: Multiple exons with same dbId for dbId "IDFMTSTR" (%p and %p)\n", Exon_getDbID(e),compEx, e);
        //exit(1);
      }
    }
*/
  }

  IDType *uniqueIds = IDHash_getKeys(exHash);
  int nUniqueId = IDHash_getNumValues(exHash);
  

  char tmpStr[1024];
  char *qStr = NULL;
  int lenNum;
  int maxSize = 16384;

  if ((qStr = (char *)calloc(655500,sizeof(char))) == NULL) {
    fprintf(stderr,"Failed allocating qStr\n");
    return NULL;
  }

  IDHash *dnaFeatIdToExHash  = IDHash_new(IDHASH_MEDIUM);
  IDHash *protFeatIdToExHash = IDHash_new(IDHASH_MEDIUM);
  IDHash *idToEx;

  for (i=0; i<nUniqueId; i+=maxSize) {
    int endPoint = sprintf(qStr, "SELECT sf.feature_type, sf.feature_id, sf.exon_id FROM supporting_feature sf WHERE sf.exon_id IN (" );

    int j;
    for (j=0; j<maxSize && j+i<nUniqueId; j++) {
      if (j!=0) {
        qStr[endPoint++] = ',';
        qStr[endPoint++] = ' ';
      }
      lenNum = sprintf(tmpStr,IDFMTSTR,uniqueIds[i+j]);
      memcpy(&(qStr[endPoint]), tmpStr, lenNum);
      endPoint+=lenNum;
    }
    qStr[endPoint++] = ')';
    qStr[endPoint] = '\0';
  
  
    StatementHandle *sth = sfa->prepare((BaseAdaptor *)sfa,qStr,strlen(qStr));
    sth->execute(sth);
  
    ResultRow *row;
    while ((row = sth->fetchRow(sth))) {
      char * type = row->getStringAt(row,0);
      IDType sfId = row->getLongLongAt(row,1);
      IDType exId = row->getLongLongAt(row,2);
  
      if (type[0] == 'd') {
        idToEx = dnaFeatIdToExHash;
      } else {
        idToEx = protFeatIdToExHash;
      }
      if (!IDHash_contains(exHash, exId)) {
        fprintf(stderr, "Warning: exon not in exhash for supporting feature [type = %s id = "IDFMTSTR" exon_id "IDFMTSTR"]\n", type, sfId, exId);
      } else {
        Vector *exVec = IDHash_getValue(idToEx, sfId);
        if (!exVec) {
          exVec = Vector_new();
          IDHash_add(idToEx, sfId, exVec);
        }
        Vector_append(exVec, IDHash_getValue(exHash, exId));
      }
    }
  
    sth->finish(sth);
  }

  free(uniqueIds);

  Vector *out = Vector_new();

  for (i=0; i<2;i++) {
    char *type;
    BaseAdaptor *ba;
    if (i==0) {
      idToEx = dnaFeatIdToExHash;
      ba = (BaseAdaptor*)DBAdaptor_getDNAAlignFeatureAdaptor(sfa->dba);
      type = "dna";
    } else {
      idToEx = protFeatIdToExHash;
      ba = (BaseAdaptor*)DBAdaptor_getProteinAlignFeatureAdaptor(sfa->dba);
      type = "protein";
    }

    if (IDHash_getNumValues(idToEx)) {
      fprintf(stderr, "Have %d ids of type %s\n", IDHash_getNumValues(idToEx), type);
  
      IDType *idArray = IDHash_getKeys(idToEx);
      Vector *idVec = Vector_new();

      int j;
      for (j=0; j < IDHash_getNumValues(idToEx); j++) {
        Vector_addElement(idVec, &idArray[j]);
      }
  
      Vector *features = BaseAdaptor_fetchAllByDbIDList((BaseAdaptor *)ba, idVec, slice);
  
      for (j=0; j<Vector_getNumElement(features); j++) {
        BaseAlignFeature *feature = Vector_getElementAt(features, j);
  
        Vector *exVec = IDHash_getValue(idToEx, BaseAlignFeature_getDbID(feature));
        //fprintf(stderr,"nexon in exVec = %d sf id = %d\n",Vector_getNumElement(exVec),BaseAlignFeature_getDbID(feature));
        int k;
        for (k=0; k<Vector_getNumElement(exVec); k++) {
          Exon *exon = Vector_getElementAt(exVec, k);
          //fprintf(stderr," exon id = "IDFMTSTR"\n", Exon_getDbID(exon));
        
          BaseAlignFeature *newFeature;
          if (slice == NULL) {
            newFeature = (BaseAlignFeature *)SeqFeature_transfer((SeqFeature *)feature, Exon_getSlice(exon));
          } else {
            newFeature = feature;
          }
  // NIY: Free feature???
          if (newFeature) {
            Exon_addSupportingFeature(exon, (SeqFeature*)newFeature);
            Vector_addElement(out, newFeature);
          } else if (slice == NULL) {
            fprintf(stderr,"Failed to transfer feature\n");
          }
        }
      }
      free(idArray);
      Vector_free(idVec);
      Vector_free(features);
    }
  }

  IDHash_free(exHash, Vector_free);
  IDHash_free(dnaFeatIdToExHash, Vector_free);
  IDHash_free(protFeatIdToExHash, Vector_free);

// Hack - if exon hasn't had any support added then it doesn't have any so set to emptyVector to flag that we've looked for
// support and didn't find any in db
  for (i=0;i<Vector_getNumElement(exons);i++) {
    Exon *exon = Vector_getElementAt(exons, i);
    if (exon->supportingFeatures == NULL) exon->supportingFeatures = emptyVector;
  }

  free(qStr);
  return out;
}


/*
=head2 store

  Arg [1]    : Int $exonsID
               The dbID of an EnsEMBL exon to associate with
               supporting features.
  Arg [2]    : Ref to array of Bio::EnsEMBL::BaseAlignFeature
               (the support)
  Example    : $sfa->store($exon_id, \@features);
  Description: Stores a set of alignment features and associates an
               EnsEMBL exon with them
  Returntype : none
  Exceptions : thrown when invalid dbID is passed to this method
  Caller     : TranscriptAdaptor
  Status     : Stable

=cut
*/
void SupportingFeatureAdaptor_store(SupportingFeatureAdaptor *sfa, IDType exonDbID, Vector *alnObjs) {
  char pepCheckSql[1024];
  char dnaCheckSql[1024];
  char assocCheckSql[1024];
  char assocWriteSql[1024];

// Note added in hcoverage so transcript_supporting_feature and supporting_feature code match
  sprintf(pepCheckSql,
      "SELECT protein_align_feature_id "
      "FROM protein_align_feature "
      "WHERE seq_region_id = %"IDFMTSTR
      " AND   seq_region_start = %%ld"
      " AND   seq_region_end   = %%ld"
      " AND   seq_region_strand = %%d"
      " AND   hit_name = '%%s'"
      " AND   hit_start = %%d"
      " AND   hit_end   = %%d"
      " AND   analysis_id = %"IDFMTSTR
      " AND   hcoverage = %%f"
      " AND   cigar_line = '%%s'");

  sprintf(dnaCheckSql,
      "SELECT dna_align_feature_id "
      "FROM dna_align_feature "
      "WHERE seq_region_id = %"IDFMTSTR
      " AND   seq_region_start = %%ld"
      " AND   seq_region_end   = %%ld"
      " AND   seq_region_strand = %%d"
      " AND   hit_name = '%%s'"
      " AND   hit_start = %%d"
      " AND   hit_end   = %%d"
      " AND   analysis_id = %"IDFMTSTR
      " AND   cigar_line = '%%s'"
      " AND   hcoverage = %%f"
      " AND   hit_strand = %%d");

  sprintf(assocCheckSql, 
      "SELECT * "
      "FROM  supporting_feature "
      "WHERE  exon_id = "IDFMTSTR
      " AND   feature_type = '%%s'"
      " AND   feature_id   = %"IDFMTSTR, exonDbID);

  sprintf(assocWriteSql, 
          "INSERT into supporting_feature (exon_id, feature_id, feature_type) values(%"IDFMTSTR", %"IDFMTSTR", '%%s')");

  StatementHandle *pepCheckSth   = sfa->prepare((BaseAdaptor *)sfa, pepCheckSql, strlen(pepCheckSql));
  StatementHandle *dnaCheckSth   = sfa->prepare((BaseAdaptor *)sfa, dnaCheckSql, strlen(dnaCheckSql));
  StatementHandle *assocCheckSth = sfa->prepare((BaseAdaptor *)sfa, assocCheckSql, strlen(assocCheckSql));
  StatementHandle *sfSth         = sfa->prepare((BaseAdaptor *)sfa, assocWriteSql, strlen(assocWriteSql));

  DNAAlignFeatureAdaptor *dnaAdaptor     = DBAdaptor_getDNAAlignFeatureAdaptor(sfa->dba);
  ProteinAlignFeatureAdaptor *pepAdaptor = DBAdaptor_getProteinAlignFeatureAdaptor(sfa->dba);
  SliceAdaptor *sliceAdaptor             = DBAdaptor_getSliceAdaptor(sfa->dba);

  int i;
  for (i=0; i < Vector_getNumElement(alnObjs); i++) {
    BaseAlignFeature *f = Vector_getElementAt(alnObjs, i);

/* I really wish I understood why this is necessary - until I do, I'm not doing these transfers
    // check that the feature is in toplevel coords
    if($f->slice->start != 1 || $f->slice->strand != 1) {
    #move feature onto a slice of the entire seq_region
      my $tls = $self->db->get_sliceAdaptor->fetch_by_region($f->slice->coord_system->name(),
                                                             $f->slice->seq_region_name(),
                                                             undef, #start
                                                             undef, #end
                                                             undef, #strand
                                                             $f->slice->coord_system->version());
      $f = $f->transfer($tls);

      if(!$f) {
        throw('Could not transfer Feature to slice of ' .
              'entire seq_region prior to storing');
      }
    }
*/

    Class_assertType(CLASS_BASEALIGNFEATURE, f->objectType);

/*
    if(!$f->isa("Bio::EnsEMBL::BaseAlignFeature")){
      throw("$f must be an align feature otherwise" .
            "it can't be stored");
    }
*/
    
    IDType sfDbID;
    char *type = NULL;
    BaseFeatureAdaptor *adap;
    StatementHandle *checkSth;

    IDType seqRegionId = SliceAdaptor_getSeqRegionId(sliceAdaptor, BaseAlignFeature_getSlice(f));
    
// Note - moved the checkSth execute into the condition because I can't do the variable args
    if (seqRegionId) {
      if (f->objectType == CLASS_DNADNAALIGNFEATURE) {
        adap     = (BaseFeatureAdaptor*)dnaAdaptor;      
        type     = "dna_align_feature";

      checkSth = dnaCheckSth;
      checkSth->execute(checkSth, seqRegionId, 
                        BaseAlignFeature_getSeqRegionStart((SeqFeature*)f), 
                        BaseAlignFeature_getSeqRegionEnd((SeqFeature*)f), 
                        BaseAlignFeature_getSeqRegionStrand((SeqFeature*)f), 
                                  BaseAlignFeature_getHitSeqName(f), 
                                  BaseAlignFeature_getHitStart(f), 
                                  BaseAlignFeature_getHitEnd(f), 
                                  Analysis_getDbID(BaseAlignFeature_getAnalysis(f)), 
                                  BaseAlignFeature_getCigarString(f), 
                                  BaseAlignFeature_gethCoverage(f), 
                                  DNAAlignFeature_getHitStrand((DNAAlignFeature *)f));

      } else if (f->objectType == CLASS_DNAPEPALIGNFEATURE) {
        adap     = (BaseFeatureAdaptor*)pepAdaptor;
        type     = "protein_align_feature";

      checkSth = pepCheckSth;
      checkSth->execute(checkSth, seqRegionId, 
                        BaseAlignFeature_getSeqRegionStart((SeqFeature*)f), 
                        BaseAlignFeature_getSeqRegionEnd((SeqFeature*)f), 
                        BaseAlignFeature_getSeqRegionStrand((SeqFeature*)f), 
                                  BaseAlignFeature_getHitSeqName(f), 
                                  BaseAlignFeature_getHitStart(f), 
                                  BaseAlignFeature_getHitEnd(f), 
                                  Analysis_getDbID(BaseAlignFeature_getAnalysis(f)), 
                                  BaseAlignFeature_getCigarString(f),
                                  BaseAlignFeature_gethCoverage(f));

      } else {
        fprintf(stderr, "Warning: Supporting feature of unknown type. Skipping\n");
        continue;
      }

      /// HOW?? - moved into conditions above
      //    checkSth->execute(@check_args);

      if (checkSth->numRows(checkSth) > 0) {
        ResultRow *row = checkSth->fetchRow(checkSth);
        sfDbID = row->getLongLongAt(row, 0);

      } else {
        Vector *vec = Vector_new();
        Vector_addElement(vec, f);
        //BaseFeatureAdaptor_store(adap, vec);
        adap->store((BaseAdaptor*)adap, vec);
        Vector_free(vec);
     
        sfDbID = BaseAlignFeature_getDbID(f);
      }

    // now check association
    assocCheckSth->execute(assocCheckSth, type, sfDbID);

    if (checkSth->numRows(assocCheckSth) == 0) {
      sfSth->execute(sfSth, exonDbID, sfDbID, type);
    }

    dnaCheckSth->finish(dnaCheckSth);
    pepCheckSth->finish(pepCheckSth);
    assocCheckSth->finish(assocCheckSth);
    sfSth->finish(sfSth);
    } else {
      fprintf(stderr, "Error getting sequence region ID for slice");
    }
  }

  return;
}
