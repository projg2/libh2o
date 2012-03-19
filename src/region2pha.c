/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include <math.h>

#include "region2.h"
#include "xmath.h"

/* Based on IF97-Rev, s. 6.3.1: The Backward Equation T(p, h) ... */

/* coefficient table; n[0] added for convenience */
static const double n[] = {
	+0.00000000000000E+00,

	+0.10898952318288E+04, +0.84951654495535E+03,
	-0.10781748091826E+03, +0.33153654801263E+02,
	-0.74232016790248E+01, +0.11765048724356E+02,
	+0.18445749355790E+01, -0.41792700549624E+01,
	+0.62478196935812E+01, -0.17344563108114E+02, /* [10] */

	-0.20058176862096E+03, +0.27196065473796E+03,
	-0.45511318285818E+03, +0.30919688604755E+04,
	+0.25226640357872E+06, -0.61707422868339E-02,
	-0.31078046629583E+00, +0.11670873077107E+02,
	+0.12812798404046E+09, -0.98554909623276E+09, /* [20] */

	+0.28224546973002E+10, -0.35948971410703E+10,
	+0.17227349913197E+10, -0.13551334240775E+05,
	+0.12848734664650E+08, +0.13865724283226E+01,
	+0.23598832556514E+06, -0.13105236545054E+08,
	+0.73999835474766E+04, -0.55196697030060E+06, /* [30] */

	+0.37154085996233E+07, +0.19127729239660E+05,
	-0.41535164835634E+06, -0.62459855192507E+02
};

static const int I[] = {
	0,

	0, 0, 0, 0, 0, 0, /* [6] */
	1, 1, 1, 1, 1, 1, 1, 1, 1, /* [15] */
	2, 2, 2, 2, 2, 2, 2, 2, /* [23] */
	3, 3,
	4, 4, 4, /* [28] */
	5, 5, 5,
	6, 6, /* [33] */
	7
};

static const double Jpows[] = {
	0, 1, 2, 3, 7, 9, 11, 12, 18, 20, 24,
	28, 32, 34, 36, 38, 40, 42, 44
};

static const int J[] = {
	0,

	0, 1, 2, 3, 4, 9, /* [6] */
	0, 1, 2, 3, 4, 5, 6, 8, 18, /* [15] */
	0, 2, 4, 14, 15, 16, 17, 18, /* [23] */
	10, 18,
	7, 12, 18, /* [28] */
	12, 14, 17,
	13, 18,
	11
};

static const double hstar = 2000; /* [kJ/kg] */

double h2o_region2a_T_ph(double p, double h) /* [MPa, kJ/kg] -> [K] */
{
	double eta = h / hstar;
	double etaexpr = eta - 2.1;

	double sum = 0;

	int i;

	double ppowers[8], etapowers[19];

	fill_powers_incr(ppowers, 8, p);
	fill_powers(etapowers, Jpows, 0, 19, etaexpr, 0);

	for (i = 1; i <= 34; ++i)
	{
		double pipow = ppowers[I[i]];
		double etapow = etapowers[J[i]];

		sum += n[i] * pipow * etapow;
	}

	return sum;
}
