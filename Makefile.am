ACLOCAL_AMFLAGS = -I m4

h2odir = $(includedir)/h2o
pkgconfigdir = $(libdir)/pkgconfig

lib_LTLIBRARIES = libh2o.la

libh2o_la_SOURCES = src/boundaries.c \
	src/region1.c src/region1ph.c src/region1ps.c \
	src/region1hs.c \
	src/region2.c src/region2sub.c src/region2meta.c \
	src/region2pha.c src/region2phb.c src/region2phc.c \
	src/region2psa.c src/region2psb.c src/region2psc.c \
	src/region3.c \
	src/region4.c src/region5.c \
	src/consts.h src/xmath.h
libh2o_la_LIBADD = -lm
libh2o_la_LDFLAGS = -version-info 0:0:0

h2o_HEADERS = src/boundaries.h \
	src/region1.h src/region2.h src/region3.h src/region4.h src/region5.h

pkgconfig_DATA = libh2o.pc

TESTS = if97-test-values
check_PROGRAMS = $(TESTS)

if97_test_values_SOURCES = tests/if97-test-values.c
if97_test_values_CPPFLAGS = -I$(top_srcdir)/src
if97_test_values_LDADD = libh2o.la
