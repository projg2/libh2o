ACLOCAL_AMFLAGS = -I m4

h2odir = $(includedir)/h2o

lib_LTLIBRARIES = libh2o.la

libh2o_la_SOURCES = src/saturation.c
libh2o_la_LIBADD = -lm
libh2o_la_LDFLAGS = -version-info 0:0:0

h2o_HEADERS = src/saturation.h

noinst_PROGRAMS = test

test_SOURCES = test.c
test_LDADD = libh2o.la

check_LTLIBRARIES = libh2o-test.la

libh2o_test_la_SOURCES = tests/tests.c
libh2o_test_la_LDFLAGS = -shared

TESTS = saturation-test
check_PROGRAMS = $(TESTS)

saturation_test_SOURCES = tests/saturation-test.c
saturation_test_CPPFLAGS = -Isrc
saturation_test_LDADD = libh2o.la libh2o-test.la
