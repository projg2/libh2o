/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif


#include "region2.h"
#include "xmath.h"

/* Based on IF97-Rev, s. 6.3.1: The Backward Equation T(p, h) ... */

/* coefficient table; n[0] added for convenience */
static const double n[] = {
	+0.00000000000000E+00,

	-0.32368398555242E+13, +0.73263350902181E+13,
	+0.35825089945447E+12, -0.58340131851590E+12,
	-0.10783068217470E+11, +0.20825544563171E+11,
	+0.61074783564516E+06, +0.85977722535580E+06,
	-0.25745723604170E+05, +0.31081088422714E+05, /* [10] */

	+0.12082315865936E+04, +0.48219755109255E+03,
	+0.37966001272486E+01, -0.10842984880077E+02,
	-0.45364172676660E-01, +0.14559115658698E-12,
	+0.11261597407230E-11, -0.17804982240686E-10,
	+0.12324579690832E-06, -0.11606921130984E-05, /* [20] */

	+0.27846367088554E-04, -0.59270038474176E-03,
	+0.12918582991878E-02
};

static const double Ipows[] = {
	-7, -6, -5, -2, -1, 0, 1, 2, 6
};

static const int I[] = {
	0,

	0, 0, 1, 1, 2, 2,
	3, 3, 4, 4, 5, 5, /* [12] */
	6, 6, 7, /* [15] */
	8, 8, 8, 8, 8, 8, 8, 8
};

static const double Jpows[] = {
	0, 1, 2, 4, 8, 10, 12, 16, 20, 22
};

static const int J[] = {
	0,

	0, 3, 0, 2, 0, 2,
	0, 1, 0, 2, 0, 1, /* [12] */
	3, 4, 3, /* [15] */
	0, 1, 3, 5, 6, 7, 8, 9
};

static const double hstar = 2000; /* [kJ/kg] */

double h2o_region2c_T_ph(double p, double h) /* [MPa, kJ/kg] -> [K] */
{
	double eta = h / hstar;
	double piexpr = p + 25;
	double etaexpr = eta - 1.8;

	double sum = 0;

	int i;

	double pipowers[9], etapowers[10];

	fill_powers(pipowers, Ipows, 5, 9, piexpr, 0);
	fill_powers(etapowers, Jpows, 0, 10, etaexpr, 0);

	for (i = 1; i <= 23; ++i)
	{
		double pipow = pipowers[I[i]];
		double etapow = etapowers[J[i]];

		sum += n[i] * pipow * etapow;
	}

	return sum;
}
