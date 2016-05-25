#
# The main EnsC makefile.am
# Process this file with automake to create Makefile.in.
#

SUBDIRS = Util DataModel DataAdaptors Compara OBDA Programs Tests

ACLOCAL_AMFLAGS = -I m4


lib_LTLIBRARIES = libensc.la

libensc_la_SOURCES = 

libensc_la_LIBADD = \
Util/libUtil.la \
DataModel/Basic/libDataModelBasic.la \
DataModel/Mapper/libDataModelMapper.la \
DataModel/libDataModel.la \
DataAdaptors/libDataAdaptors.la \
Compara/DataModel/libComparaDataModel.la \
Compara/DataAdaptors/libComparaDataAdaptors.la \
OBDA/libOBDA.la \
$(NULL)


# keep libtool up to date automatically.
LIBTOOL_DEPS = @LIBTOOL_DEPS@
libtool: $(LIBTOOL_DEPS)
	$(SHELL) ./config.status libtool


# Extra files to remove for the maintainer-clean target.
# Note you cannot use this target to remove directories,
# hence the extra "local" target.
#
MAINTAINERCLEANFILES = \
$(top_srcdir)/Makefile.in \
$(top_srcdir)/aclocal.m4 \
$(top_srcdir)/config.h.in \
$(top_srcdir)/config.h.in~ \
$(top_srcdir)/configure \
$(top_srcdir)/ensc_version.m4

GENERATED_FILES = $(top_srcdir)/config $(top_srcdir)/autom4te.cache

maintainer-clean-local:
	rm -rf $(GENERATED_FILES)
	mkdir $(top_srcdir)/config
	rm -rf $(top_srcdir)/build/*
	git checkout $(top_srcdir)/build
