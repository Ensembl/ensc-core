# EnsC makefile

CC = cc

MYSQLCFLAGS=${shell mysql_config --cflags}
MYSQLLIBS=${shell mysql_config --libs}

CFLAGS =  -g -IDataModel -IDataAdaptors -IUtil -IDataModel/Mapper -IDataModel/Basic -ICompara/DataAdaptors -ICompara/DataModel $(MYSQLCFLAGS)
LDFLAGS = $(MYSQLLIBS)

BINDIR=../bin

all	: directories \
          libensc.a \
          Programs/testdbc \
          Programs/countsnps \
          Tests/CacheTest \
          Tests/CigarStrUtilTest \
          Tests/ClassTest \
          Tests/ComparaDBAdaptorTest \
          Tests/ComparaDNAAlignFeatureAdaptorTest \
          Tests/DBAdaptorTest \
          Tests/DNAAlignFeatureTest \
          Tests/DNAPepAlignFeatureTest \
          Tests/EcoStringTest \
          Tests/HomologyTest \
          Tests/RepeatFeatureTest \
          Tests/SeqUtilTest \
          Tests/SimpleFeatureTest \
          Tests/SliceAdaptorTest \
          Tests/StreamTest \
          Tests/StrUtilTest \
          Tests/SyntenyTest \
          Tests/TranslateTest \
          Tests/VectorTest

directories	:
	mkdir -p objectfiles
	mkdir -p objectfiles/DataAdaptors
	mkdir -p objectfiles/DataModel
	mkdir -p objectfiles/Util
	mkdir -p objectfiles/Programs
	mkdir -p objectfiles/Tests
	mkdir -p objectfiles/DataModel/Basic
	mkdir -p objectfiles/DataModel/Mapper
	mkdir -p objectfiles/Compara/DataModel
	mkdir -p objectfiles/Compara/DataAdaptors

objectfiles/%.o: %.c
	$(CC) -c $(CFLAGS) -o $@ $<

libensc.a	: \
	objectfiles/Compara/DataAdaptors/BaseComparaAdaptor.o \
	objectfiles/Compara/DataAdaptors/ComparaDBAdaptor.o \
	objectfiles/Compara/DataAdaptors/ComparaDNAAlignFeatureAdaptor.o \
	objectfiles/Compara/DataAdaptors/DNAFragAdaptor.o \
	objectfiles/Compara/DataAdaptors/GenomeDBAdaptor.o \
	objectfiles/Compara/DataAdaptors/GenomicAlignAdaptor.o \
	objectfiles/Compara/DataAdaptors/HomologyAdaptor.o \
	objectfiles/Compara/DataAdaptors/SyntenyAdaptor.o \
	objectfiles/Compara/DataModel/DNAFrag.o \
	objectfiles/Compara/DataModel/GenomeDB.o \
	objectfiles/Compara/DataModel/GenomicAlign.o \
	objectfiles/Compara/DataModel/Homology.o \
	objectfiles/Compara/DataModel/SyntenyRegion.o \
	objectfiles/DataAdaptors/AnalysisAdaptor.o \
	objectfiles/DataAdaptors/AssemblyMapperAdaptor.o \
	objectfiles/DataAdaptors/BaseAdaptor.o \
	objectfiles/DataAdaptors/BaseFeatureAdaptor.o \
	objectfiles/DataAdaptors/ChromosomeAdaptor.o \
	objectfiles/DataAdaptors/CloneAdaptor.o \
	objectfiles/DataAdaptors/DBAdaptor.o \
	objectfiles/DataAdaptors/DBConnection.o \
	objectfiles/DataAdaptors/DBEntryAdaptor.o \
	objectfiles/DataAdaptors/DNAAlignFeatureAdaptor.o \
	objectfiles/DataAdaptors/ExonAdaptor.o \
	objectfiles/DataAdaptors/GeneAdaptor.o \
	objectfiles/DataAdaptors/MetaContainer.o \
	objectfiles/DataAdaptors/MysqlResultRow.o \
	objectfiles/DataAdaptors/MysqlStatementHandle.o \
	objectfiles/DataAdaptors/PredictionTranscriptAdaptor.o \
	objectfiles/DataAdaptors/ProteinAlignFeatureAdaptor.o \
	objectfiles/DataAdaptors/RawContigAdaptor.o \
	objectfiles/DataAdaptors/RepeatConsensusAdaptor.o \
	objectfiles/DataAdaptors/RepeatFeatureAdaptor.o \
	objectfiles/DataAdaptors/SequenceAdaptor.o \
	objectfiles/DataAdaptors/SimpleFeatureAdaptor.o \
	objectfiles/DataAdaptors/SliceAdaptor.o \
	objectfiles/DataAdaptors/SupportingFeatureAdaptor.o \
	objectfiles/DataAdaptors/TranscriptAdaptor.o \
	objectfiles/DataAdaptors/TranslationAdaptor.o \
	objectfiles/DataModel/Analysis.o \
	objectfiles/DataModel/AnnotatedSeqFeature.o \
	objectfiles/DataModel/AssemblyMapper.o \
	objectfiles/DataModel/BaseAlignFeature.o \
	objectfiles/DataModel/BaseContig.o \
	objectfiles/DataModel/Basic/Vector.o \
	objectfiles/DataModel/Chromosome.o \
	objectfiles/DataModel/CigarStrUtil.o \
	objectfiles/DataModel/Class.o \
	objectfiles/DataModel/Clone.o \
	objectfiles/DataModel/DBEntry.o \
	objectfiles/DataModel/DNAAlignFeature.o \
	objectfiles/DataModel/DNAPepAlignFeature.o \
	objectfiles/DataModel/EnsC.o \
	objectfiles/DataModel/Exon.o \
	objectfiles/DataModel/FeaturePair.o \
	objectfiles/DataModel/FeatureSet.o \
	objectfiles/DataModel/Gene.o \
	objectfiles/DataModel/GenomicRange.o \
	objectfiles/DataModel/IdentityXref.o \
	objectfiles/DataModel/Mapper.o \
	objectfiles/DataModel/Mapper/MapperCoordinate.o \
	objectfiles/DataModel/Mapper/MapperGap.o \
	objectfiles/DataModel/Mapper/MapperPair.o \
	objectfiles/DataModel/Mapper/MapperPairSet.o \
	objectfiles/DataModel/Mapper/MapperRangeSet.o \
	objectfiles/DataModel/Mapper/MapperUnit.o \
	objectfiles/DataModel/PredictionTranscript.o \
	objectfiles/DataModel/RawContig.o \
	objectfiles/DataModel/RepeatConsensus.o \
	objectfiles/DataModel/RepeatFeature.o \
	objectfiles/DataModel/SeqFeature.o \
	objectfiles/DataModel/SeqFeatureFactory.o \
	objectfiles/DataModel/Sequence.o \
	objectfiles/DataModel/SimpleFeature.o \
	objectfiles/DataModel/Slice.o \
	objectfiles/DataModel/Species.o \
	objectfiles/DataModel/StableIdInfo.o \
	objectfiles/DataModel/StickyExon.o \
	objectfiles/DataModel/Transcript.o \
	objectfiles/DataModel/Translation.o \
	objectfiles/Util/CHash.o \
	objectfiles/Util/Cache.o \
	objectfiles/Util/EcoString.o \
	objectfiles/Util/Error.o \
	objectfiles/Util/FileUtil.o \
	objectfiles/Util/IDHash.o \
	objectfiles/Util/MysqlUtil.o \
	objectfiles/Util/ProcUtil.o \
	objectfiles/Util/SeqUtil.o \
	objectfiles/Util/StrUtil.o \
	objectfiles/Util/Stream.o \
	objectfiles/Util/StringHash.o \
	objectfiles/Util/translate.o
	ar rcv libensc.a objectfiles/DataModel/*.o objectfiles/DataAdaptors/*.o objectfiles/Util/*.o \
               objectfiles/DataModel/Basic/*.o objectfiles/DataModel/Mapper/*.o \
               objectfiles/Compara/DataAdaptors/*.o objectfiles/Compara/DataModel/*.o
	ranlib libensc.a

Programs/testdbc: \
	libensc.a \
	Programs/testdbc.c
	$(CC) $(CFLAGS) -o Programs/testdbc Programs/testdbc.c -L. -lensc $(LDFLAGS)

Programs/countsnps: \
	libensc.a \
	Programs/countsnps.c
	$(CC) $(CFLAGS) -o Programs/countsnps Programs/countsnps.c -L. -lensc $(LDFLAGS)

Tests/CacheTest: \
	libensc.a \
	Tests/CacheTest.c
	$(CC) $(CFLAGS) -o Tests/CacheTest Tests/CacheTest.c -L. -lensc $(LDFLAGS)

Tests/ClassTest: \
	libensc.a \
	Tests/ClassTest.c
	$(CC) $(CFLAGS) -o Tests/ClassTest Tests/ClassTest.c -L. -lensc $(LDFLAGS)

Tests/CigarStrUtilTest: \
	libensc.a \
	Tests/CigarStrUtilTest.c
	$(CC) $(CFLAGS) -o Tests/CigarStrUtilTest Tests/CigarStrUtilTest.c -L. -lensc $(LDFLAGS)

Tests/ComparaDBAdaptorTest: \
	libensc.a \
	Tests/ComparaDBAdaptorTest.c
	$(CC) $(CFLAGS) -o Tests/ComparaDBAdaptorTest Tests/ComparaDBAdaptorTest.c -L. -lensc $(LDFLAGS)

Tests/ComparaDNAAlignFeatureAdaptorTest: \
	libensc.a \
	Tests/ComparaDNAAlignFeatureAdaptorTest.c
	$(CC) $(CFLAGS) -o Tests/ComparaDNAAlignFeatureAdaptorTest Tests/ComparaDNAAlignFeatureAdaptorTest.c -L. -lensc $(LDFLAGS)

Tests/DBAdaptorTest: \
	libensc.a \
	Tests/DBAdaptorTest.c
	$(CC) $(CFLAGS) -o Tests/DBAdaptorTest Tests/DBAdaptorTest.c -L. -lensc $(LDFLAGS)

Tests/DNAAlignFeatureTest: \
	libensc.a \
	Tests/DNAAlignFeatureTest.c
	$(CC) $(CFLAGS) -o Tests/DNAAlignFeatureTest Tests/DNAAlignFeatureTest.c -L. -lensc $(LDFLAGS)

Tests/DNAPepAlignFeatureTest: \
	libensc.a \
	Tests/DNAPepAlignFeatureTest.c
	$(CC) $(CFLAGS) -o Tests/DNAPepAlignFeatureTest Tests/DNAPepAlignFeatureTest.c -L. -lensc $(LDFLAGS)

Tests/EcoStringTest: \
	libensc.a \
	Tests/EcoStringTest.c
	$(CC) $(CFLAGS) -o Tests/EcoStringTest Tests/EcoStringTest.c -L. -lensc $(LDFLAGS)

Tests/HomologyTest: \
	libensc.a \
	Tests/HomologyTest.c
	$(CC) $(CFLAGS) -o Tests/HomologyTest Tests/HomologyTest.c -L. -lensc $(LDFLAGS)

Tests/RepeatFeatureTest: \
	libensc.a \
	Tests/RepeatFeatureTest.c
	$(CC) $(CFLAGS) -o Tests/RepeatFeatureTest Tests/RepeatFeatureTest.c -L. -lensc $(LDFLAGS)

Tests/SeqUtilTest: \
	libensc.a \
	Tests/SeqUtilTest.c
	$(CC) $(CFLAGS) -o Tests/SeqUtilTest Tests/SeqUtilTest.c -L. -lensc $(LDFLAGS)

Tests/SimpleFeatureTest: \
	libensc.a \
	Tests/SimpleFeatureTest.c
	$(CC) $(CFLAGS) -o Tests/SimpleFeatureTest Tests/SimpleFeatureTest.c -L. -lensc $(LDFLAGS)

Tests/SliceAdaptorTest: \
	libensc.a \
	Tests/SliceAdaptorTest.c
	$(CC) $(CFLAGS) -o Tests/SliceAdaptorTest Tests/SliceAdaptorTest.c -L. -lensc $(LDFLAGS)

Tests/StreamTest: \
	libensc.a \
	Tests/StreamTest.c
	$(CC) $(CFLAGS) -o Tests/StreamTest Tests/StreamTest.c -L. -lensc $(LDFLAGS)

Tests/StrUtilTest: \
	libensc.a \
	Tests/StrUtilTest.c
	$(CC) $(CFLAGS) -o Tests/StrUtilTest Tests/StrUtilTest.c -L. -lensc $(LDFLAGS)

Tests/SyntenyTest: \
	libensc.a \
	Tests/SyntenyTest.c
	$(CC) $(CFLAGS) -o Tests/SyntenyTest Tests/SyntenyTest.c -L. -lensc $(LDFLAGS)

Tests/TranslateTest: \
	libensc.a \
	Tests/TranslateTest.c
	$(CC) $(CFLAGS) -o Tests/TranslateTest Tests/TranslateTest.c -L. -lensc $(LDFLAGS)

Tests/VectorTest: \
	libensc.a \
	Tests/VectorTest.c
	$(CC) $(CFLAGS) -o Tests/VectorTest Tests/VectorTest.c -L. -lensc $(LDFLAGS)

# deletes compiled class files
clean	:	
	rm -f libensc.a
	rm -f Programs/testdbc
	rm -f Programs/countsnps
	find ./ -name "*.o" -exec rm {} \;
	find ./Tests ! -name "*.c" ! -name "*.h" -type f -perm -0700 -exec rm {} \;

#DO NOT REMOVE THIS LINE OR ADD ANYTHING AFTER IT - ADDMAKE DEPENDS ON IT FOR MAKING DEPENDENCIES

