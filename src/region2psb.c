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

/* Based on IF97-Rev, s. 6.3.2: The Backward Equation T(p, s) ... */

/* coefficient table; n[0] added for convenience */
static const double n[] = {
	+0.00000000000000E+00,

	+0.31687665083497E+6, +0.20864175881858E+2,
	-0.39859399803599E+6, -0.21816058518877E+2,
	+0.22369785194242E+6, -0.27841703445817E+4,
	+0.99207436071480E+1, -0.75197512299157E+5,
	+0.29708605951158E+4, -0.34406878548526E+1, /* [10] */

	+0.38815564249115E+0, +0.17511295085750E+5,
	-0.14237112854449E+4, +0.10943803364167E+1,
	+0.89971619308495E+0, -0.33759740098958E+4,
	+0.47162885818355E+3, -0.19188241993679E+1,
	+0.41078580492196E+0, -0.33465378172097E+0, /* [20] */

	+0.13870034777505E+4, -0.40663326195838E+3,
	+0.41727347159610E+2, +0.21932549434532E+1,
	-0.10320050009077E+1, +0.35882943516703E+0,
	+0.52511453726066E-2, +0.12838916450705E+2,
	-0.28642437219381E+1, +0.56912683664855E+0, /* [30] */

	-0.99962954584931E-1, -0.32632037778459E-2,
	+0.23320922576723E-3, -0.15334809857450E+0,
	+0.29072288239902E-1, +0.37534702741167E-3,
	+0.17296691702411E-2, -0.38556050844504E-3,
	-0.35017712292608E-4, -0.14566393631492E-4, /* [40] */

	+0.56420857267269E-5, +0.41286150074605E-7,
	-0.20684671118824E-7, +0.16409393674725E-8
};

static const int I[] = {
	0,

	-6, -6, -5, -5, -4, -4, -4,
	-3, -3, -3, -3,
	-2, -2, -2, -2,
	-1, -1, -1, -1, -1,
	0, 0, 0, 0, 0, 0, 0,
	1, 1, 1, 1, 1, 1,
	2, 2, 2, 3, 3, 3,
	4, 4, 5, 5, 5
};

static const int J[] = {
	0,

	0, 11, 0, 11, 0, 1, 11,
	0, 1, 11, 12,
	0, 1, 6, 10,
	0, 1, 5, 8, 9,
	0, 1, 2, 4, 5, 6, 9,
	0, 1, 2, 3, 7, 8,
	0, 1, 5, 0, 1, 3,
	0, 1, 0, 1, 2
};

static const double sstar = 0.7853; /* [kJ/kgK] */

double h2o_region2b_T_ps(double p, double s) /* [MPa, kJ/kgK] -> [K] */
{
	double sigma = s / sstar;
	double sigmaexpr = 10 - sigma;

	double sum = 0;

	int i;

	double pipowers_store[6+6], sigmapowers[13];
	double* pipowers = &pipowers_store[6];

	fill_powers_incr(pipowers, 6, p);
	fill_powers_decr(pipowers, -6, p);
	fill_powers_incr(sigmapowers, 13, sigmaexpr);

	for (i = 1; i <= 44; ++i)
	{
		double pipow = pipowers[I[i]];
		double sigmapow = sigmapowers[J[i]];

		sum += n[i] * pipow * sigmapow;
	}

	return sum;
}
