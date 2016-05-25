
include $(top_srcdir)/Makefile.am.common

SUBDIRS = .

noinst_LTLIBRARIES = libComparaDataAdaptors.la

include_HEADERS = \
BaseComparaAdaptor.h \
ComparaAdaptorTypes.h \
ComparaDBAdaptor.h \
ComparaDNAAlignFeatureAdaptor.h \
DNAFragAdaptor.h \
GenomeDBAdaptor.h \
GenomicAlignAdaptor.h \
HomologyAdaptor.h \
SyntenyAdaptor.h \
$(NULL)

libComparaDataAdaptors_la_SOURCES = \
BaseComparaAdaptor.c \
ComparaDBAdaptor.c \
ComparaDNAAlignFeatureAdaptor.c \
DNAFragAdaptor.c \
GenomeDBAdaptor.c \
GenomicAlignAdaptor.c \
HomologyAdaptor.c \
SyntenyAdaptor.c \
$(NULL)


libComparaDataAdaptors_la_LIBADD = \
$(top_builddir)/Util/libUtil.la \
$(top_builddir)/DataModel/Basic/libDataModelBasic.la \
$(top_builddir)/DataModel/Mapper/libDataModelMapper.la \
$(top_builddir)/DataModel/libDataModel.la \
$(top_builddir)/DataAdaptors/libDataAdaptors.la \
$(top_builddir)/Compara/DataModel/libComparaDataModel.la \
$(NULL)


# Extra files to remove for the maintainer-clean target.
#
MAINTAINERCLEANFILES = $(top_srcdir)/Compara/DataAdaptors/Makefile.in

