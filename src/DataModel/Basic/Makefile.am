
include $(top_srcdir)/Makefile.am.common

SUBDIRS = .

noinst_LTLIBRARIES = libDataModelBasic.la

pkginclude_HEADERS = \
Vector.h \
IDHash.h \
$(NULL)

libDataModelBasic_la_SOURCES = \
Vector.c \
IDHash.c \
$(NULL)

libDataModelBasic_la_LIBADD = \
$(top_builddir)/Util/libUtil.la \
$(NULL)


# Extra files to remove for the maintainer-clean target.
#
MAINTAINERCLEANFILES = $(top_srcdir)/DataModel/Basic/Makefile.in

