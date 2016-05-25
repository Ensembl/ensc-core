
include $(top_srcdir)/Makefile.am.common

SUBDIRS = .

noinst_LTLIBRARIES = libComparaDataModel.la

include_HEADERS = \
ComparaDataModelTypes.h \
DNAFrag.h \
GenomeDB.h \
GenomicAlign.h \
Homology.h \
SyntenyRegion.h \
$(NULL)

libComparaDataModel_la_SOURCES = \
DNAFrag.c \
GenomeDB.c \
GenomicAlign.c \
Homology.c \
SyntenyRegion.c \
$(NULL)


libComparaDataModel_la_LIBADD = \
$(top_builddir)/Util/libUtil.la \
$(top_builddir)/DataModel/Basic/libDataModelBasic.la \
$(top_builddir)/DataModel/Mapper/libDataModelMapper.la \
$(top_builddir)/DataModel/libDataModel.la \
$(top_builddir)/DataAdaptors/libDataAdaptors.la \
$(NULL)


# Extra files to remove for the maintainer-clean target.
#
MAINTAINERCLEANFILES = $(top_srcdir)/Compara/DataModel/Makefile.in

