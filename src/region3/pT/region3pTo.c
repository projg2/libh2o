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

	+0.128746023979718E-34, -0.735234770382342E-11,
	+0.289078692149150E-02, +0.244482731907223E+00,
	+0.141733492030985E-23, -0.354533853059476E-28,
	-0.594539202901431E-17, -0.585188401782779E-08,

	+0.201377325411803E-05, +0.138647388209306E+01,
	-0.173959365084772E-04, +0.137680878349369E-02,
	+0.814897605805513E-14, +0.425596631351839E-25,
	-0.387449113787755E-17, +0.139814747930240E-12,

	-0.171849638951521E-02, +0.641890529513296E-21,
	+0.118960578072018E-10, -0.155282762571611E-17,
	+0.233907907347507E-07, -0.174093247766213E-12,
	+0.377682649089149E-08, -0.516720236575302E-10
};

static const double Ipows[] = {
	0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 14, 20, 24
};

static const int I[] = {
	0,

	0, 0, 0, 2, 3, 4, 4, 4, 4, 4, 5, 5, 6, 7, 8, 8, 8, 9, 9, 10, 10, 11, 11, 12
};

static const double Jpows[] = {
	-12, -10, -8, -5, -4, -3, -1, 0, 1
};

static const int J[] = {
	0,

	0, 4, 6, 6, 1, 0, 2, 3, 4, 6, 4, 5, 2, 0, 1, 2, 4, 0, 2, 0, 2, 0, 1, 0
};

static const double vstar = 0.0034; /* [m³/kg] */
static const double pstar = 23; /* [MPa] */
static const double Tstar = 650; /* [K] */

double h2o_region3o_v_pT(double p, double T)
{
	double pi = p / pstar;
	double theta = T / Tstar;

	return twoarg_poly_value(sqrt(pi - 0.974), theta - 0.996,
			I, Ipows, 0, 13, 0,
			J, Jpows, 7, 9, 0,
			n, 24) * vstar;
}
