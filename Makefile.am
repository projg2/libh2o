ACLOCAL_AMFLAGS = -I m4

h2odir = $(includedir)/h2o

lib_LTLIBRARIES = libh2o.la

libh2o_la_SOURCES = src/lowlevel/saturation.c
libh2o_la_LIBADD = -lm
libh2o_la_LDFLAGS = -version-info 0:0:0

h2o_HEADERS = src/lowlevel/saturation.h

noinst_PROGRAMS = test

test_SOURCES = src/test.c
test_LDADD = libh2o.la
