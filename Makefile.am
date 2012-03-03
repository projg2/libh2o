ACLOCAL_AMFLAGS = -I m4

lib_LTLIBRARIES = libh2o.la

libh2o_la_SOURCES = src/lowlevel/saturation.c
libh2o_la_LIBADD = -lm
libh2o_la_LDFLAGS = -version-info 0:0:0

noinst_PROGRAMS = test

test_SOURCES = src/test.c
test_LDADD = libh2o.la
