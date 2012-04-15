/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include "region3.h"
#include "xmath.h"

/* Supplementary Release on Backward Equations for Specific Volume
 * as a Function of Pressure and Temperature v(p,T)
 * for Region 3 of the IAPWS Industrial Formulation 1997 for the
 * Thermodynamic Properties of Water and Steam */

/* coefficient table; n[0] added for convenience */
static const double n[] = {
	+0.000000000000000E+00,

	+0.260702058647537E+10, -0.188277213604704E+15,
	+0.554923870289667E+19, -0.758966946387758E+23,
	+0.413865186848908E+27, -0.815038000738060E+12,
	-0.381458260489955E+33, -0.123239564600519E-01,

	+0.226095631437174E+08, -0.495017809506720E+12,
	+0.529482996422863E+16, -0.444359478746295E+23,
	+0.521635864527315E+35, -0.487095672740742E+55,
	-0.714430209937547E+06, +0.127868634615495E+00,

	-0.100752127917598E+02, +0.777451437960990E+07,
	-0.108105480796471E+25, -0.357578581169659E-05,
	-0.212857169423484E+01, +0.270706111085238E+30,
	-0.695953622348829E+33, +0.110609027472280E+00,

	+0.721559163361354E+02, -0.306367307532219E+15,
	+0.265839618885530E-04, +0.253392392889754E-01,
	-0.214443041836579E+03, +0.937846601489667E+00,
	+0.223184043101700E+01, +0.338401222509191E+02,

	+0.494237237179718E+21, -0.198068404154428E+00,
	-0.141415349881140E+31, -0.993862421613651E+02,
	+0.125070534142731E+03, -0.996473529004439E+03,
	+0.473137909872765E+05, +0.116662121219322E+33,

	-0.315874976271533E+16, -0.445703369196945E+33,
	+0.642794932373694E+33
};

static const double Ipows[] = {
	-12, -10, -8, -6, -5, -4, -3, -2, -1, 0, 1, 2, 4, 5, 6, 10, 14
};

static const int I[] = {
	0,

	0, 0, 0, 0, 0, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 4, 4, 5, 5, 6, 6, 6, 6,
	7, 7, 7, 8, 8, 8, 9, 9, 9, 9, 10, 10, 11, 12, 13, 13, 14, 15, 15, 16
};

static const double Jpows[] = {
	0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16, 18, 20, 22, 24, 36
};

static const int J[] = {
	0,

	11, 12, 13, 14, 15, 11, 16, 6, 9, 10, 11, 13, 16, 17, 8, 4, 5, 7,
	12, 1, 3, 13, 14, 2, 3, 9, 0, 1, 3, 0, 1, 2, 10, 0, 12, 1, 0, 0, 1, 11,
	4, 10, 9
};

static const double vstar = 0.0026; /* [m³/kg] */
static const double pstar = 24; /* [MPa] */
static const double Tstar = 650; /* [K] */

double h2o_region3l_v_pT(double p, double T)
{
	double pi = p / pstar;
	double theta = T / Tstar;

	double sum = poly_value(pi - 0.908, theta - 0.989,
			I, Ipows, 9, 17, 0,
			J, Jpows, 0, 18, 0,
			n, 43);

	return pow4(sum) * vstar;
}
