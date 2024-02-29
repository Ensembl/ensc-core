/*
 * Copyright [1999-2024] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::DBSQL::ProteinAlignFeatureAdaptor - 
Adaptor for ProteinAlignFeatures
*/
#include "ProteinAlignFeatureAdaptor.h"

#include "AnalysisAdaptor.h"
#include "DNAPepAlignFeature.h"
#include "SliceAdaptor.h"
#include "ChainedAssemblyMapper.h"



NameTableType ProteinAlignFeatureAdaptor_tableNames = {{"protein_align_feature","paf"},
                                                       {"external_db", "exdb"},
                                                       {NULL, NULL}};

NameTableType ProteinAlignFeatureAdaptor_leftJoins = {{"external_db","exdb.external_db_id = paf.external_db_id"},
                                                      {NULL,NULL}};


ProteinAlignFeatureAdaptor *ProteinAlignFeatureAdaptor_new(DBAdaptor *dba) {
  ProteinAlignFeatureAdaptor *pafa;

  if ((pafa = (ProteinAlignFeatureAdaptor *)calloc(1,sizeof(ProteinAlignFeatureAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for ProteinAlignFeatureAdaptor\n");
    return NULL;
  }
  BaseFeatureAdaptor_init((BaseFeatureAdaptor *)pafa, dba, PROTEINALIGNFEATURE_ADAPTOR);

  pafa->getTables                  = ProteinAlignFeatureAdaptor_getTables;
  pafa->getColumns                 = ProteinAlignFeatureAdaptor_getColumns;
  pafa->store                      = (BaseAdaptor_StoreFunc)ProteinAlignFeatureAdaptor_store;
  pafa->objectsFromStatementHandle = (BaseAdaptor_ObjectsFromStatementHandleFunc)ProteinAlignFeatureAdaptor_objectsFromStatementHandle;
  pafa->leftJoin                   = ProteinAlignFeatureAdaptor_leftJoin;

  return pafa;
}

NameTableType *ProteinAlignFeatureAdaptor_getTables(void) {
  return &ProteinAlignFeatureAdaptor_tableNames;
}



char *ProteinAlign_cols[] = {
             "paf.protein_align_feature_id",
             "paf.seq_region_id",
             "paf.seq_region_start",
             "paf.seq_region_end",
             "paf.analysis_id",
             "paf.seq_region_strand",
             "paf.hit_start",
             "paf.hit_end",
             "paf.hit_name",
             "paf.cigar_line",
             "paf.evalue",
             "paf.perc_ident",
             "paf.score",
             "paf.external_db_id",
             "paf.hcoverage",
	     "exdb.db_name",
	     "exdb.db_display_name",
             NULL };

char **ProteinAlignFeatureAdaptor_getColumns(void) {
  return ProteinAlign_cols;
}


NameTableType *ProteinAlignFeatureAdaptor_leftJoin(void) {
  return &ProteinAlignFeatureAdaptor_leftJoins;
}


/*
=head2 store

  Arg [1]    : list of Bio::EnsEMBL::DnaPepAlignFeature @feats
  Example    : $protein_align_feature_adaptor->store(@feats);
  Description: stores a list of ProteinAlignFeatures in the database
  Returntype : none
  Exceptions : throw if any of the provided features cannot be stored
               which may occur if:
                 * The feature does not have an associated Slice
                 * The feature does not have an associated analysis
                 * The Slice the feature is associated with is on a seq_region
                   unknown to this database
              A warning is given if:
                 * The feature has already been stored in this db
  Caller     : Pipeline
  Status     : Stable

=cut
*/

int ProteinAlignFeatureAdaptor_store(BaseFeatureAdaptor *bfa, Vector *features) {
  if (features == NULL || Vector_getNumElement(features) == 0) {
    fprintf(stderr,"Must call store with features\n");
    exit(1);
  }

  NameTableType *tables = bfa->getTables();
  char *tableName = (*tables)[0][NAME];

  DBAdaptor *db                    = bfa->dba;
  AnalysisAdaptor *analysisAdaptor = DBAdaptor_getAnalysisAdaptor(db);
  SliceAdaptor *sliceAdaptor = DBAdaptor_getSliceAdaptor(db);

  char qStr[1024];
  sprintf(qStr, "INSERT INTO %s (seq_region_id, seq_region_start, seq_region_end,"
                                "seq_region_strand, hit_start, hit_end,"
                                "hit_name, cigar_line,"
                                "analysis_id, score, evalue, perc_ident, external_db_id, hcoverage) "
                       "VALUES (%" IDFMTSTR ",%%d,%%d,%%d,%%d,%%d,\'%%s\',\'%%s\',%"IDFMTSTR ",%%f,%%f,%%f,%" IDFMTSTR ",%%f)", tableName);

  StatementHandle *sth = bfa->prepare((BaseAdaptor *)bfa, qStr,strlen(qStr));

  int i;
  for (i=0; i<Vector_getNumElement(features); i++) {
    char fixedCigar[1024];

    DNAPepAlignFeature *feat = Vector_getElementAt(features, i);

    if (feat == NULL) {
      fprintf(stderr, "feature is NULL in ProteinAlignFeature_store\n");
      exit(1);
    }

    Class_assertType(CLASS_DNAPEPALIGNFEATURE, feat->objectType);

    if (DNAPepAlignFeature_isStored(feat, db)) {
      fprintf(stderr, "DNAPepAlignFeature ["IDFMTSTR"] is already stored in this database.\n", DNAPepAlignFeature_getDbID(feat) );
      continue;
    }

    BaseFeatureAdaptor_checkStartEndStrand(bfa,
                                           DNAPepAlignFeature_getHitStart(feat),
                                           DNAPepAlignFeature_getHitEnd(feat),
                                           1,
                                           NULL);



    char *cigarString = DNAPepAlignFeature_getCigarString(feat);

    if (cigarString == NULL) {
      sprintf(fixedCigar, "%ldM", DNAPepAlignFeature_getLength(feat));
      cigarString = fixedCigar;
      fprintf(stderr, "DNAPepAlignFeature does not define a cigar_string.\n"
                      "Assuming ungapped block with cigar_line = %s.\n", cigarString);
    }

    if (DNAPepAlignFeature_getHitSeqName(feat) == NULL) {
      fprintf(stderr, "DNAPepAlignFeature must define an hseqname.\n");
      exit(1);
    }

    Analysis *analysis = DNAPepAlignFeature_getAnalysis(feat);
    if (analysis == NULL) {
      fprintf(stderr,"An analysis must be attached to the features to be stored.\n");
      exit(1);
    }

    // store the analysis if it has not been stored yet
    if (Analysis_isStored(analysis, db)) {
      AnalysisAdaptor_store(analysisAdaptor, analysis);
    }

/*
   my $original = $feat;
   my $seq_region_id;
   ($feat, $seq_region_id) = $self->_pre_store($feat);
*/
    IDType seqRegionId = BaseFeatureAdaptor_preStore(bfa, (SeqFeature*)feat);
// Note using SeqRegionStart etc here rather than Start - should have same effect as perl's transfer
    sth->execute(sth, (IDType)seqRegionId,
                 DNAPepAlignFeature_getSeqRegionStart((SeqFeature*)feat),
                 DNAPepAlignFeature_getSeqRegionEnd((SeqFeature*)feat),
                 DNAPepAlignFeature_getSeqRegionStrand((SeqFeature*)feat),
                      DNAPepAlignFeature_getHitStart(feat),
                      DNAPepAlignFeature_getHitEnd(feat),
                      DNAPepAlignFeature_getHitSeqName(feat),
                      cigarString,
                      (IDType)Analysis_getDbID(analysis),
                      DNAPepAlignFeature_getScore(feat),
                      DNAPepAlignFeature_getpValue(feat),
                      DNAPepAlignFeature_getPercId(feat),
                      (IDType)DNAPepAlignFeature_getExternalDbID(feat),
                      DNAPepAlignFeature_gethCoverage(feat));


    DNAPepAlignFeature_setDbID(feat,sth->getInsertId(sth));
    if (!DNAPepAlignFeature_getDbID(feat)) {
      // Not sure it's correct to finish here, but it's probably better than not doing it
      sth->finish(sth);
      exit(1);
    }
    DNAPepAlignFeature_setAdaptor(feat, (BaseAdaptor *)bfa);
  }

  sth->finish(sth);

  return 1;
}


/*
=head2 _objs_from_sth

  Arg [1]    : DBI statement handle $sth
               an exectuted DBI statement handle generated by selecting 
               the columns specified by _columns() from the table specified 
               by _table()
  Example    : @dna_dna_align_feats = $self->_obj_from_hashref
  Description: PROTECTED implementation of superclass abstract method. 
               Creates DnaDnaAlignFeature objects from a DBI hashref
  Returntype : listref of Bio::EnsEMBL::ProteinAlignFeatures
  Exceptions : none
  Caller     : Bio::EnsEMBL::BaseFeatureAdaptor::generic_fetch
  Status     : Stable

=cut
*/
Vector *ProteinAlignFeatureAdaptor_objectsFromStatementHandle(BaseFeatureAdaptor *bfa,
                                                           StatementHandle *sth,
                                                           AssemblyMapper *assMapper,
                                                           Slice *destSlice) {
  SliceAdaptor *sa     = DBAdaptor_getSliceAdaptor(bfa->dba);
  AnalysisAdaptor *aa  = DBAdaptor_getAnalysisAdaptor(bfa->dba);

  Vector *features = Vector_new();
  IDHash *sliceHash = IDHash_new(IDHASH_SMALL);

  long         destSliceStart;
  long         destSliceEnd;
  int          destSliceStrand;
  long         destSliceLength;
  char *       destSliceSrName;
  IDType       destSliceSrId = 0;

  if (destSlice) {
    destSliceStart  = Slice_getStart(destSlice);
    destSliceEnd    = Slice_getEnd(destSlice);
    destSliceStrand = Slice_getStrand(destSlice);
    destSliceLength = Slice_getLength(destSlice);
    destSliceSrName = Slice_getSeqRegionName(destSlice);
    destSliceSrId   = Slice_getSeqRegionId(destSlice);
  }

  ResultRow *row;
  while ((row = sth->fetchRow(sth))) {
    IDType proteinAlignFeatureId = row->getLongLongAt(row,0);
    IDType seqRegionId           = row->getLongLongAt(row,1);
    long seqRegionStart          = row->getLongAt(row,2);
    long seqRegionEnd            = row->getLongAt(row,3);
    IDType analysisId            = row->getLongLongAt(row,4);
    int seqRegionStrand          = row->getIntAt(row,5);
    int hitStart                 = row->getIntAt(row,6);
    int hitEnd                   = row->getIntAt(row,7);
    char *hitName                = row->getStringAt(row,8);
    char *cigarLine              = row->getStringAt(row,9);
    double eValue                = row->getDoubleAt(row,10);
    double percIdent             = row->getDoubleAt(row,11);
    double score                 = row->getDoubleAt(row,12);
    IDType externalDbId          = row->getLongLongAt(row,13);
    double hCoverage             = row->getDoubleAt(row,14);
    char *externalDbName         = row->getStringAt(row,15);
    char *externalDbDisplayName  = row->getStringAt(row,16);

    // get the analysis object
    Analysis *analysis = AnalysisAdaptor_fetchByDbID(aa, analysisId);

    if (! IDHash_contains(sliceHash, seqRegionId)) {
      IDHash_add(sliceHash, seqRegionId, SliceAdaptor_fetchBySeqRegionId(sa, seqRegionId, POS_UNDEF, POS_UNDEF, STRAND_UNDEF));
    }
    Slice *slice = IDHash_getValue(sliceHash, seqRegionId);

    Slice *pafSlice = slice;

    char *srName      = Slice_getSeqRegionName(slice);
    CoordSystem *srCs = Slice_getCoordSystem(slice);

    //
    // remap the feature coordinates to another coord system
    // if a mapper was provided
    //
    if (assMapper != NULL) {
      MapperRangeSet *mrs;

      // Slightly suspicious about need for this if statement so left in perl statements for now
      if (destSlice != NULL &&
          assMapper->objectType == CLASS_CHAINEDASSEMBLYMAPPER) {
        mrs = ChainedAssemblyMapper_map(assMapper, srName, seqRegionStart, seqRegionEnd, seqRegionStrand, srCs, 1, destSlice);
      } else {
        mrs = AssemblyMapper_fastMap(assMapper, srName, seqRegionStart, seqRegionEnd, seqRegionStrand, srCs, NULL);
      }

      // skip features that map to gaps or coord system boundaries
      if (MapperRangeSet_getNumRange(mrs) == 0) {
        continue;
      }
      MapperRange *range = MapperRangeSet_getRangeAt(mrs, 0);
      if (range->rangeType == MAPPERRANGE_GAP) {
        fprintf(stderr,"Got a mapper gap in gene obj_from_sth - not sure if this is allowed\n");
        exit(1);
      } else {
        MapperCoordinate *mc = (MapperCoordinate *)range;

        seqRegionId     = mc->id;
        seqRegionStart  = mc->start;
        seqRegionEnd    = mc->end;
        seqRegionStrand = mc->strand;
      }

      MapperRangeSet_free(mrs);

      if (! IDHash_contains(sliceHash, seqRegionId)) {
        IDHash_add(sliceHash, seqRegionId, SliceAdaptor_fetchBySeqRegionId(sa, seqRegionId, POS_UNDEF, POS_UNDEF, STRAND_UNDEF));
      }
      pafSlice = IDHash_getValue(sliceHash, seqRegionId);
    }

    if (destSlice != NULL) {
      if (destSliceStart != 1 || destSliceStrand != 1) {
        if (destSliceStrand == 1) {
          seqRegionStart = seqRegionStart - destSliceStart + 1;
          seqRegionEnd   = seqRegionEnd - destSliceStart + 1;
        } else {
          long tmpSeqRegionStart = seqRegionStart;
          seqRegionStart = destSliceEnd - seqRegionEnd + 1;
          seqRegionEnd   = destSliceEnd - tmpSeqRegionStart + 1;

          seqRegionStrand = -seqRegionStrand;
        }
      }
      // throw away features off the end of the requested slice
      if (seqRegionEnd < 1 || seqRegionStart > destSliceLength || (destSliceSrId != seqRegionId)) {
        continue;
      }
      pafSlice = destSlice;
    }

    // Finally, create the new DnaAlignFeature.
    DNAPepAlignFeature *paf = DNAPepAlignFeature_new();

    DNAPepAlignFeature_setSlice(paf,pafSlice);
    DNAPepAlignFeature_setStart(paf, seqRegionStart);
    DNAPepAlignFeature_setEnd(paf, seqRegionEnd);
    DNAPepAlignFeature_setStrand(paf, seqRegionStrand);

    DNAPepAlignFeature_setHitSeqName(paf,hitName);
    DNAPepAlignFeature_setHitStart(paf, hitStart);
    DNAPepAlignFeature_setHitEnd(paf, hitEnd);
    DNAPepAlignFeature_setHitStrand(paf, 1);

    DNAPepAlignFeature_setScore(paf, score);
    DNAPepAlignFeature_setpValue(paf, eValue);
    DNAPepAlignFeature_setPercId(paf, percIdent);

    DNAPepAlignFeature_setCigarString(paf, cigarLine);

    DNAPepAlignFeature_setAnalysis(paf,analysis);
    DNAPepAlignFeature_setAdaptor(paf, (BaseAdaptor *)bfa);
    DNAPepAlignFeature_setDbID(paf, proteinAlignFeatureId);

    DNAPepAlignFeature_setExternalDbID(paf, externalDbId);
    DNAPepAlignFeature_sethCoverage(paf, hCoverage);
    DNAPepAlignFeature_setDbName(paf, externalDbName);
    DNAPepAlignFeature_setDbDisplayName(paf, externalDbDisplayName);


    Vector_addElement(features, paf);
  }

  IDHash_free(sliceHash, NULL);
  return features;
}


/*
=head2 list_dbIDs

  Arg [1]    : none
  Example    : @feature_ids = @{$protein_align_feature_adaptor->list_dbIDs()};
  Description: Gets an array of internal ids for all protein align 
               features in the current db
  Arg[1]     : <optional> int. not 0 for the ids to be sorted by the seq_region.
  Returntype : listref of ints
  Exceptions : none
  Caller     : ?
  Status     : Stable

=cut
*/
Vector *ProteinAlignFeatureAdaptor_listDbIDs(ProteinAlignFeatureAdaptor *pafa, int ordered) {

// NIY: Shouldn't really do this direct BaseAdaptor call but ...
  return BaseAdaptor_listDbIDs((BaseAdaptor *)pafa, "protein_align_feature", NULL, ordered);
}

