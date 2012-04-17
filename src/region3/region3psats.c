/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include "region3.h"
#include "xmath.h"

/* Revised Supplementary Release on Backward Equations for the Functions
 * T(p,h), v(p,h) and T(p,s), v(p,s) for Region 3 of the IAPWS Industrial
 * Formulation 1997 for the Thermodynamic Properties of Water and Steam
 * s. 4.3: Boundary Equations psat(h) and psat(s) */

/* coefficient table; n[0] added for convenience */
static const double n[] = {
	+0.000000000000000E+00,

	+0.639767553612785E00, -0.129727445396014E02,
	-0.224595125848403E16, +0.177466741801846E07,
	+0.717079349571538E10, -0.378829107169011E18,
	-0.955586736431328E35, +0.187269814676188E24,

	+0.119254746466473E12, +0.110649277244882E37
};

static const double Ipows[] = {
	0, 1, 4, 12, 16, 24, 28, 32
};

static const int I[] = {
	0,

	0, 1, 1, 2, 3, 3, 4, 5, 6, 7
};

static const double Jpows[] = {
	0, 1, 4, 7, 10, 14, 18, 32, 36
};

static const int J[] = {
	0,

	0, 1, 7, 3, 2, 5, 8, 4, 0, 6
};

static const double pstar = 22; /* [MPa] */
static const double sstar = 5.2; /* [kJ/kgK] */

double h2o_region3_psat_s(double s)
{
	double sigma = s / sstar;

	return twoarg_poly_value(sigma - 1.03, sigma - 0.699,
			I, Ipows, 0, 8, 0,
			J, Jpows, 0, 9, 0,
			n, 10) * pstar;
}
