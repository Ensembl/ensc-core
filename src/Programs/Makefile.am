
include $(top_srcdir)/Makefile.am.common

SUBDIRS = .

bin_PROGRAMS = 

if HAVE_SAMTOOLS
  bin_PROGRAMS += bamcount bamcount_exon bamcov bammap countsnps fastasplit indicate moveperl testdbc

  bamcount_SOURCES = bamcount.c
  bamcount_exon_SOURCES = bamcount_exon.c
  bamcov_SOURCES = bamcov.c
  bammap_SOURCES = bammap.c
  countsnps_SOURCES = countsnps.c
  fastasplit_SOURCES = fastasplit.c
  indicate_SOURCES = indicate.c indicate.h
  moveperl_SOURCES = moveperl.c
  testdbc_SOURCES = testdbc.c

  if HAVE_LIBTCMALLOC
    bin_PROGRAMS += RefineSolexaGenes

    RefineSolexaGenes_SOURCES = RefineSolexaGenes.c RefineSolexaGenes.h 
  endif
endif



# Extra files to remove for the maintainer-clean target.
#
MAINTAINERCLEANFILES = $(top_srcdir)/Programs/Makefile.in

