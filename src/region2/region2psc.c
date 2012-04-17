/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif


#include "region2.h"
#include "xmath.h"

/* Based on IF97-Rev, s. 6.3.2: The Backward Equation T(p, s) ... */

/* coefficient table; n[0] added for convenience */
static const double n[] = {
	+0.00000000000000E+00,

	+0.90968501005365E+03, +0.24045667088420E+04,
	-0.59162326387130E+03, +0.54145404128074E+03,
	-0.27098308411192E+03, +0.97976525097926E+03,
	-0.46966772959435E+03, +0.14399274604723E+02,
	-0.19104204230429E+02, +0.53299167111971E+01, /* [10] */

	-0.21252975375934E+02, -0.31147334413760E+00,
	+0.60334840894623E+00, -0.42764839702509E-01,
	+0.58185597255259E-02, -0.14597008284753E-01,
	+0.56631175631027E-02, -0.76155864584577E-04,
	+0.22440342919332E-03, -0.12561095013413E-04, /* [20] */

	+0.63323132660934E-06, -0.20541989675375E-05,
	+0.36405370390082E-07, -0.29759897789215E-08,
	+0.10136618529763E-07, +0.59925719692351E-11,
	-0.20677870105164E-10, -0.20874278181886E-10,
	+0.10162166825089E-09, -0.16429828281347E-09
};

static const double Ipows[] = {
	-2, -1, 0, 1, 2, 3, 4, 5, 6, 7
};

static const int I[] = {
	0,

	0, 0,
	1,
	2, 2, 2, 2,
	3, 3, 3, 3,
	4, 4, 4,
	5, 5, 5,
	6, 6, 6,
	7, 7, 7,
	8, 8,
	9, 9, 9, 9, 9
};

static const double Jpows[] = {
	0, 1, 2, 3, 4, 5
};

static const int J[] = {
	0,

	0, 1,
	0,
	0, 1, 2, 3,
	0, 1, 3, 4,
	0, 1, 2,
	0, 1, 5,
	0, 1, 4,
	0, 1, 2,
	0, 1,
	0, 1, 3, 4, 5
};

static const double sstar = 2.9251; /* [kJ/kgK] */

double h2o_region2c_T_ps(double p, double s) /* [MPa, kJ/kgK] -> [K] */
{
	double sigma = s / sstar;

	return twoarg_poly_value(p, 2 - sigma,
			I, Ipows, 2, 10, 0,
			J, Jpows, 0, 6, 0,
			n, 30);
}
