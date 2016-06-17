
include $(top_srcdir)/Makefile.am.common

SUBDIRS = Basic Mapper .

noinst_LTLIBRARIES = libDataModel.la

include_HEADERS = \
Analysis.h \
AnnotatedSeqFeature.h \
AssemblyMapper.h \
Attribute.h \
BaseAlignFeature.h \
BaseAssemblyMapper.h \
BaseContig.h \
ChainedAssemblyMapper.h \
Chromosome.h \
CigarStrUtil.h \
Class.h \
Clone.h \
CoordSystem.h \
DBEntry.h \
DNAAlignFeature.h \
DNAPepAlignFeature.h \
DataModelTypes.h \
EnsRoot.h \
Exon.h \
FeaturePair.h \
FeatureSet.h \
Gene.h \
GenomicRange.h \
IdentityXref.h \
Intron.h \
IntronSupportingEvidence.h \
Mapper.h \
MapperDefs.h \
Object.h \
PredictionExon.h \
PredictionTranscript.h \
ProjectionSegment.h \
RawContig.h \
RepeatConsensus.h \
RepeatFeature.h \
SeqEdit.h \
SeqFeature.h \
SeqFeatureFactory.h \
SeqRegionRange.h \
Sequence.h \
SimpleFeature.h \
Slice.h \
Species.h \
StableIdInfo.h \
StickyExon.h \
Storable.h \
TopLevelAssemblyMapper.h \
Transcript.h \
Translation.h \
$(NULL)

libDataModel_la_SOURCES = \
Analysis.c \
AnnotatedSeqFeature.c \
AssemblyMapper.c \
Attribute.c \
BaseAlignFeature.c \
BaseAssemblyMapper.c \
BaseContig.c \
ChainedAssemblyMapper.c \
Chromosome.c \
CigarStrUtil.c \
Class.c \
Clone.c \
CoordSystem.c \
DBEntry.c \
DNAAlignFeature.c \
DNAPepAlignFeature.c \
Exon.c \
FeaturePair.c \
FeatureSet.c \
Gene.c \
GenomicRange.c \
IdentityXref.c \
Intron.c \
IntronSupportingEvidence.c \
Mapper.c \
Object.c \
PredictionExon.c \
PredictionTranscript.c \
ProjectionSegment.c \
RawContig.c \
RepeatConsensus.c \
RepeatFeature.c \
SeqEdit.c \
SeqFeature.c \
SeqFeatureFactory.c \
SeqRegionRange.c \
Sequence.c \
SimpleFeature.c \
Slice.c \
Species.c \
StableIdInfo.c \
Storable.c \
TopLevelAssemblyMapper.c \
Transcript.c \
Translation.c \
$(NULL)


libDataModel_la_LIBADD = \
$(top_builddir)/Util/libUtil.la \
$(top_builddir)/DataModel/Basic/libDataModelBasic.la \
$(top_builddir)/DataModel/Mapper/libDataModelMapper.la \
$(NULL)


# Extra files to remove for the maintainer-clean target.
#
MAINTAINERCLEANFILES = $(top_srcdir)/DataModel/Makefile.in

