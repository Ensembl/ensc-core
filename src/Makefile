# EnsC makefile

CC = cc

MYSQLCFLAGS=${shell mysql_config --cflags}
MYSQLLIBS=${shell mysql_config --libs}

CFLAGS = -O2 -IDataModel -IDataAdaptors -IUtil -IDataModel/Mapper -IDataModel/Basic $(MYSQLCFLAGS)
LDFLAGS = $(MYSQLLIBS)

BINDIR=../bin

all	: directories \
          libensc.a \
          Programs/testdbc \
          Programs/countsnps \
          Tests/CacheTest \
          Tests/ClassTest \
          Tests/DBAdaptorTest \
          Tests/EcoStringTest \
          Tests/SeqUtilTest \
          Tests/StreamTest \
          Tests/StrUtilTest

directories	:
	mkdir -p objectfiles
	mkdir -p objectfiles/DataAdaptors
	mkdir -p objectfiles/DataModel
	mkdir -p objectfiles/Util
	mkdir -p objectfiles/Programs
	mkdir -p objectfiles/Tests
	mkdir -p objectfiles/DataModel/Basic
	mkdir -p objectfiles/DataModel/Mapper

objectfiles/%.o: %.c
	$(CC) -c $(CFLAGS) -o $@ $<

libensc.a	: \
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
	objectfiles/DataModel/Class.o \
	objectfiles/DataModel/Clone.o \
	objectfiles/DataModel/DBEntry.o \
	objectfiles/DataModel/DNAAlignFeature.o \
	objectfiles/DataModel/DNAPepAlignFeature.o \
	objectfiles/DataModel/EnsC.o \
	objectfiles/DataModel/Exon.o \
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
	objectfiles/DataModel/SimpleFeature.o \
	objectfiles/DataModel/Slice.o \
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
	objectfiles/Util/SeqUtil.o \
	objectfiles/Util/StrUtil.o \
	objectfiles/Util/Stream.o \
	objectfiles/Util/StringHash.o
	ar rcv libensc.a objectfiles/DataModel/*.o objectfiles/DataAdaptors/*.o objectfiles/Util/*.o \
               objectfiles/DataModel/Basic/*.o objectfiles/DataModel/Mapper/*.o
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

Tests/DBAdaptorTest: \
	libensc.a \
	Tests/DBAdaptorTest.c
	$(CC) $(CFLAGS) -o Tests/DBAdaptorTest Tests/DBAdaptorTest.c -L. -lensc $(LDFLAGS)

Tests/EcoStringTest: \
	libensc.a \
	Tests/EcoStringTest.c
	$(CC) $(CFLAGS) -o Tests/EcoStringTest Tests/EcoStringTest.c -L. -lensc $(LDFLAGS)

Tests/SeqUtilTest: \
	libensc.a \
	Tests/SeqUtilTest.c
	$(CC) $(CFLAGS) -o Tests/SeqUtilTest Tests/SeqUtilTest.c -L. -lensc $(LDFLAGS)

Tests/StreamTest: \
	libensc.a \
	Tests/StreamTest.c
	$(CC) $(CFLAGS) -o Tests/StreamTest Tests/StreamTest.c -L. -lensc $(LDFLAGS)

Tests/StrUtilTest: \
	libensc.a \
	Tests/StrUtilTest.c
	$(CC) $(CFLAGS) -o Tests/StrUtilTest Tests/StrUtilTest.c -L. -lensc $(LDFLAGS)

# deletes compiled class files
clean	:	
	rm -f libensc.a
	rm -f Programs/testdbc
	rm -f Programs/countsnps
	find ./ -name "*.o" -exec rm {} \;
	find ./Tests ! -name "*.c" ! -name "*.h" -type f -perm -0700 -exec rm {} \;

#DO NOT REMOVE THIS LINE OR ADD ANYTHING AFTER IT - ADDMAKE DEPENDS ON IT FOR MAKING DEPENDENCIES

