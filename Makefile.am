ACLOCAL_AMFLAGS = -I m4

h2odir = $(includedir)/h2o

lib_LTLIBRARIES = libh2o.la

libh2o_la_SOURCES = src/saturation.c src/region1.c src/region1ph.c
libh2o_la_LIBADD = -lm
libh2o_la_LDFLAGS = -version-info 0:0:0

h2o_HEADERS = src/saturation.h src/region1.h

check_LTLIBRARIES = libh2o-test.la

libh2o_test_la_SOURCES = tests/tests.c
libh2o_test_la_LDFLAGS = -shared

TESTS = if97-test-values
check_PROGRAMS = $(TESTS)

if97_test_values_SOURCES = tests/if97-test-values.c
if97_test_values_CPPFLAGS = -Isrc
if97_test_values_LDADD = libh2o.la libh2o-test.la
