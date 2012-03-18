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

static const int I[] = {
	0,

	-7, -7, -6, -6, -5, -5,
	-2, -2, -1, -1, 0, 0, /* [12] */
	1, 1, 2, /* [15] */
	6, 6, 6, 6, 6, 6, 6, 6
};

static const int J[] = {
	0,

	0, 4, 0, 2, 0, 2,
	0, 1, 0, 2, 0, 1, /* [12] */
	4, 8, 4, /* [15] */
	0, 1, 4, 10, 12, 16, 20, 22
};

static const double hstar = 2000; /* [kJ/kg] */

double h2o_region2c_T_ph(double p, double h) /* [MPa, kJ/kg] -> [K] */
{
	double eta = h / hstar;
	double piexpr = p + 25;
	double etaexpr = eta - 1.8;

	double sum = 0;

	int i;

	double pipowers_store[7+7], etapowers[3];
	double* pipowers = &pipowers_store[7];

	pipowers[0] = 1;
	pipowers[1] = piexpr;
	pipowers[2] = pipowers[1] * piexpr;
	pipowers[6] = pow(piexpr, 6);

	pipowers[-1] = 1 / piexpr;
	pipowers[-2] = pipowers[-1] / piexpr;
	pipowers[-5] = pow(piexpr, -5);
	pipowers[-6] = pipowers[-5] / piexpr;
	pipowers[-7] = pipowers[-6] / piexpr;

	etapowers[0] = 1;
	etapowers[1] = etaexpr;
	etapowers[2] = etaexpr * etaexpr;
	etapowers[4] = etapowers[2] * etapowers[2];

	for (i = 1; i <= 23; ++i)
	{
		double pipow = pipowers[I[i]];
		double etapow = J[i] <= 4 ? etapowers[J[i]]
			: pow(etaexpr, J[i]);

		sum += n[i] * pipow * etapow;
	}

	return sum;
}
