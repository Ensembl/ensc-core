
include $(top_srcdir)/Makefile.am.common

SUBDIRS = . 

noinst_LTLIBRARIES = libDataModelMapper.la

include_HEADERS = \
CoordPair.h \
IndelCoordinate.h \
MapperCoordinate.h \
MapperGap.h \
MapperPair.h \
MapperPairSet.h \
MapperRange.h \
MapperRangeSet.h \
MapperUnit.h \
RangeRegistry.h \
$(NULL)

libDataModelMapper_la_SOURCES = \
CoordPair.c \
IndelCoordinate.c \
MapperCoordinate.c \
MapperGap.c \
MapperPair.c \
MapperPairSet.c \
MapperRangeSet.c \
MapperUnit.c \
RangeRegistry.c \
$(NULL)


# Extra files to remove for the maintainer-clean target.
#
MAINTAINERCLEANFILES = $(top_srcdir)/DataModel/Mapper/Makefile.in

