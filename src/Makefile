# EnsC makefile

CC = cc

MYSQLCFLAGS=${shell mysql_config --cflags}
MYSQLLIBS=${shell mysql_config --libs}

CFLAGS = -IDataModel -IDataAdaptors -IUtil -IDataModel/Mapper -IDataModel/Basic $(MYSQLCFLAGS)

BINDIR=../bin

all	: libensc.a Tests/EcoStringTest

objectfiles/%.o: %.c
	$(CC) -c $(CFLAGS) -o $@ $<

libensc.a	: \
	objectfiles/DataAdaptors/AnalysisAdaptor.o \
	objectfiles/DataAdaptors/AssemblyMapperAdaptor.o \
	objectfiles/DataAdaptors/BaseAdaptor.o \
	objectfiles/DataAdaptors/ChromosomeAdaptor.o \
	objectfiles/DataAdaptors/CloneAdaptor.o \
	objectfiles/DataAdaptors/DBAdaptor.o \
	objectfiles/DataAdaptors/DBConnection.o \
	objectfiles/DataAdaptors/ExonAdaptor.o \
	objectfiles/DataAdaptors/GeneAdaptor.o \
	objectfiles/DataAdaptors/MetaContainer.o \
	objectfiles/DataAdaptors/RawContigAdaptor.o \
	objectfiles/DataAdaptors/SequenceAdaptor.o \
	objectfiles/DataAdaptors/SliceAdaptor.o \
	objectfiles/DataAdaptors/StatementHandle.o \
	objectfiles/DataAdaptors/TranscriptAdaptor.o \
	objectfiles/DataAdaptors/TranslationAdaptor.o \
	objectfiles/DataModel/Analysis.o \
	objectfiles/DataModel/AssemblyMapper.o \
	objectfiles/DataModel/BaseAlignFeature.o \
	objectfiles/DataModel/Basic/Set.o \
	objectfiles/DataModel/Chromosome.o \
	objectfiles/DataModel/Clone.o \
	objectfiles/DataModel/Exon.o \
	objectfiles/DataModel/FeatureSet.o \
	objectfiles/DataModel/Gene.o \
	objectfiles/DataModel/GenomicRange.o \
	objectfiles/DataModel/Mapper.o \
	objectfiles/DataModel/Mapper/MapperCoordinate.o \
	objectfiles/DataModel/Mapper/MapperGap.o \
	objectfiles/DataModel/Mapper/MapperPair.o \
	objectfiles/DataModel/Mapper/MapperPairSet.o \
	objectfiles/DataModel/Mapper/MapperRangeSet.o \
	objectfiles/DataModel/Mapper/MapperUnit.o \
	objectfiles/DataModel/RawContig.o \
	objectfiles/DataModel/SeqFeature.o \
	objectfiles/DataModel/Slice.o \
	objectfiles/DataModel/StableIdInfo.o \
	objectfiles/DataModel/StickyExon.o \
	objectfiles/DataModel/Transcript.o \
	objectfiles/DataModel/Translation.o \
	objectfiles/Programs/testdbc.o \
	objectfiles/Util/CHash.o \
	objectfiles/Util/EcoString.o \
	objectfiles/Util/Error.o \
	objectfiles/Util/FileUtil.o \
	objectfiles/Util/IDHash.o \
	objectfiles/Util/MysqlUtil.o \
	objectfiles/Util/SeqUtil.o \
	objectfiles/Util/StrUtil.o \
	objectfiles/Util/Stream.o \
	objectfiles/Util/StringHash.o
	ar rcv libensc.a objectfiles/DataModel/*.o objectfiles/DataAdaptors/*.o objectfiles/Util/*.o
	ranlib libensc.a

Tests/EcoStringTest: \
	libensc.a \
	Tests/EcoStringTest.c
	cc $(CFLAGS) -o Tests/EcoStringTest Tests/EcoStringTest.c -L. -lensc

# deletes compiled class files
clean	:	
	find ./ -name "*.o" -exec rm {} \;

#DO NOT REMOVE THIS LINE OR ADD ANYTHING AFTER IT - ADDMAKE DEPENDS ON IT FOR MAKING DEPENDENCIES

