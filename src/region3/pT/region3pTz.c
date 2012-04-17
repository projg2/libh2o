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

	+0.244007892290650E-10, -0.463057430331242E+07,
	+0.728803274777712E+10, +0.327776302858856E+16,
	-0.110598170118409E+10, -0.323899915729957E+13,
	+0.923814007023245E+16, +0.842250080413712E-12,

	+0.663221436245506E+12, -0.167170186672139E+15,
	+0.253749358701391E+04, -0.819731559610523E-20,
	+0.328380587890663E+12, -0.625004791171543E+08,
	+0.803197957462023E+21, -0.204397011338353E-10,

	-0.378391047055938E+04, +0.972876545938620E-02,
	+0.154355721681459E+02, -0.373962862928643E+04,
	-0.682859011374572E+11, -0.248488015614543E-03,
	+0.394536049497068E+07
};

static const double Ipows[] = {
	-8, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 6, 8
};

static const int I[] = {
	0,

	0, 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 6, 7, 8, 9, 10, 10, 11, 11, 11, 11, 12, 12
};

static const double Jpows[] = {
	-8, -6, -5, -4, -2, -1, 0, 1, 2, 3, 5, 6, 8
};

static const int J[] = {
	0,

	9, 11, 11, 12, 10, 11, 12, 4, 10, 11, 8, 1, 9, 7, 11, 1, 4, 1, 2, 3, 5, 0, 3
};

static const double vstar = 0.0038; /* [m³/kg] */
static const double pstar = 22; /* [MPa] */
static const double Tstar = 650; /* [K] */

double h2o_region3z_v_pT(double p, double T)
{
	double pi = p / pstar;
	double theta = T / Tstar;

	double sum = twoarg_poly_value(pi - 0.993, theta - 0.994,
			I, Ipows, 7, 13, 0,
			J, Jpows, 6, 13, 0,
			n, 23);

	return pow4(sum) * vstar;
}
