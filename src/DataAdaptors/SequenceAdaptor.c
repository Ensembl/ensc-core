#include "SequenceAdaptor.h"
#include "BaseAdaptor.h"
#include "MysqlUtil.h"
#include "SeqUtil.h"
#include "RawContig.h"
#include "AssemblyMapperAdaptor.h"
#include "RawContigAdaptor.h"
#include "BaseContig.h"
#include "Slice.h"
#include "StrUtil.h"

#include "StatementHandle.h"
#include "ResultRow.h"


SequenceAdaptor *SequenceAdaptor_new(DBAdaptor *dba) {
  SequenceAdaptor *sa;

  if ((sa = (SequenceAdaptor *)calloc(1,sizeof(SequenceAdaptor))) == NULL) {
    fprintf(stderr, "ERROR: Failed allocating space for SequenceAdaptor\n");
    return NULL;
  }
  BaseAdaptor_init((BaseAdaptor *)sa, dba, SEQUENCE_ADAPTOR);

  return sa;
}

char *SequenceAdaptor_fetchByRawContigStartEndStrand(SequenceAdaptor *sa, 
                                                   RawContig *rc,
                                                   int start,
                                                   int end,
                                                   char strand) {
  char qStr[256];
  StatementHandle *sth;
  ResultRow *row;

  if( start < 1 ) {
    fprintf(stderr,"ERROR: Wrong parameters to SequenceAdaptor fetch\n");
    return NULL;
  }

  if( end == -1 ) {

    sprintf(qStr,
            "SELECT c.length, SUBSTRING( d.sequence, %d )"
              " FROM dna d, contig c"
              " WHERE d.dna_id = c.dna_id"
              "  AND c.contig_id = %d", start, RawContig_getDbID(rc));

  } else {
    int length = end - start + 1;
    if( length < 1 ) {
      fprintf(stderr,"ERROR: Wrong parameters to SequenceAdaptor fetch\n");
      return NULL;
    }

    sprintf(qStr,"SELECT c.length,"
                 " SUBSTRING( d.sequence, %d, %d )"
                 " FROM dna d, contig c"
                 " WHERE d.dna_id = c.dna_id"
                 " AND c.contig_id = %d",start,length,RawContig_getDbID(rc));
  }

  if(  DBAdaptor_getDNADBAdaptor(sa->dba) ) {
    sth = DBAdaptor_prepare(DBAdaptor_getDNADBAdaptor(sa->dba),qStr,
                                                          strlen(qStr));
  } else {
    sth = sa->prepare( (BaseAdaptor *)sa, qStr, strlen(qStr) );
  }

  sth->execute(sth);

  row = sth->fetchRow(sth);

  if( row != NULL) {
    int length   = row->getIntAt(row,0);
    char *seqStr = row->getStringAt(row,1);

    /* Is this necessary ????? $seq =~ s/\s//g; */
    sth->finish(sth);
    if( strand == -1 ) {
      return SeqUtil_reverseComplement( seqStr, strlen(seqStr) );
    } else {
      return seqStr;
    }
  } else {
    sth->finish(sth);
    return NULL;
  }
}

char *SequenceAdaptor_fetchBySliceStartEndStrand(SequenceAdaptor *sa,
                                                 Slice *slice, int start, int end,
                                                 int strand) {
  char *seq;

  if (!slice ) {
    fprintf(stderr,"ERROR: need a slice to work\n");
    exit(1);
  }

  if (BaseContig_getContigType(slice) != SLICE) {
    fprintf(stderr,"ERROR: slice fetch called with something that isn't a slice\n");
    exit(1);
  }

  if (end == -1) {
    end = Slice_getChrEnd(slice) - Slice_getChrStart(slice) + 1;
  }

  // need to check the strand'edness of the slice as this
  // affects the direction in which the dna seq is grabbed
  if (Slice_getStrand(slice) == 1) {
    seq = SequenceAdaptor_fetchByAssemblyLocation(sa,
            Slice_getChrStart(slice)+start-1,
            Slice_getChrStart(slice)+end-1,
            strand,
            Slice_getChrId(slice),
            Slice_getAssemblyType(slice)
           );
  } else if (Slice_getStrand(slice) == -1 ) {
    seq = SequenceAdaptor_fetchByAssemblyLocation(sa,
            Slice_getChrEnd(slice)-end+1,
            Slice_getChrEnd(slice)-start+1,
            strand * -1, //have to make strand relative to slice's strand
            Slice_getChrId(slice),
            Slice_getAssemblyType(slice)
           );
  } else {
    fprintf(stderr,"ERROR: Incorrect strand set on slice\n");
    exit(1);
  }
  return seq;
}

char *SequenceAdaptor_fetchByAssemblyLocation(SequenceAdaptor *sa,
          int chrStart, int chrEnd, int strand, int64 chrId, char *assemblyType) {

  AssemblyMapperAdaptor *ama = DBAdaptor_getAssemblyMapperAdaptor(sa->dba);
  RawContigAdaptor *ra = DBAdaptor_getRawContigAdaptor(sa->dba);
  AssemblyMapper *assMapper = AssemblyMapperAdaptor_fetchByType(ama, assemblyType);
  char *seq;
  MapperRangeSet *coordSet;
  int i;


  coordSet = AssemblyMapper_mapCoordinatesToRawContig(assMapper, chrId, 
                                                      chrStart, chrEnd, strand );

  // for each of the pieces get sequence
  seq = StrUtil_copyString(&seq, "", 0);
  for (i=0; i<MapperRangeSet_getNumRange(coordSet); i++) {
    MapperRange *segment = MapperRangeSet_getRangeAt(coordSet,i);

    if (segment->rangeType == MAPPERRANGE_COORD) {
      MapperCoordinate *mc = (MapperCoordinate *)segment;
      RawContig *contig = RawContigAdaptor_fetchByDbID(ra, mc->id);
      char *contigSeq = SequenceAdaptor_fetchByRawContigStartEndStrand(sa, contig,
                                      mc->start, mc->end, mc->strand); 
      seq = StrUtil_appendString(seq,contigSeq);
      free(contigSeq);
    } else {
      // its a gap
      int length = segment->end - segment->start + 1;
      seq = SeqUtil_addNs(seq,length);
    }
  }

  return seq;
}

