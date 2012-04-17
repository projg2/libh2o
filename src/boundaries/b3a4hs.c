/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include "boundaries.h"
#include "xmath.h"

/* Supplementary Release on Backward Equations p(h,s) for Region 3, Equations as
 * a Function of h and s for the Region Boundaries, and an Equation Tsat(h,s)
 * for Region 4 of the IAPWS Industrial Formulation 1997 for the Thermodynamic
 * Properties of Water and Steam
 * s. 4.3: Equations h'1(s) and h'3a(s) for the Saturated Liquid Line */

/* coefficient table; n[0] added for convenience */
static const double n[] = {
	+0.00000000000000E+0,

	+0.822673364673336E+0, +0.181977213534479E+0,
	-0.112000260313624E-1, -0.746778287048033E-3,
	-0.179046263257381E+0, +0.424220110836657E-1,
	-0.341355823438768E+0, -0.209881740853565E+1,

	-0.822477343323596E+1, -0.499684082076008E+1,
	+0.191413958471069E+0, +0.581062241093136E-1,
	-0.165505498701029E+4, +0.158870443421201E+4,
	-0.850623535172818E+2, -0.317714386511207E+5,

	-0.945890406632871E+5, -0.139273847088690E-5,
	+0.631052532240980E+0
};

static const double Ipows[] = {
	0, 1, 2, 3, 4, 5, 6, 7, 10, 32
};

static const int I[] = {
	0,

	0, 0, 0, 0, 2, 3, 4, 4, 5, 5,
	6, 7, 7, 7, 8, 8, 8, 9, 9
};

static const double Jpows[] = {
	0, 1, 2, 3, 4, 6, 10, 14, 16, 20, 28, 32, 36
};

static const int J[] = {
	0,

	1, 4, 6, 8, 1, 12, 3, 8, 9, 12,
	4, 2, 10, 11, 7, 11, 12, 0, 5
};

static const double hstar = 1700; /* [kJ/kg] */
static const double sstar = 3.8; /* [kJ/kgK] */

double h2o_b3a4_h_s(double s)
{
	double sigma = s / sstar;

	return twoarg_poly_value(sigma - 1.09, sigma + 0.366E-4,
			I, Ipows, 0, 10, 0,
			J, Jpows, 0, 13, 0,
			n, 19) * hstar;
}
