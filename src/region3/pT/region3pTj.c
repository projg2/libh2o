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

	-0.111371317395540E-03, +0.100342892423685E+01,
	+0.530615581928979E+01, +0.179058760078792E-05,
	-0.728541958464774E-03, -0.187576133371704E+02,
	+0.199060874071849E-02, +0.243574755377290E+02,

	-0.177040785499444E-03, -0.259680385227130E-02,
	-0.198704578406823E+03, +0.738627790224287E-04,
	-0.236264692844138E-02, -0.161023121314333E+01,
	+0.622322971786473E+04, -0.960754116701669E-08,

	-0.510572269720488E-10, +0.767373781404211E-02,
	+0.663855469485254E-14, -0.717590735526745E-09,
	+0.146564542926508E-04, +0.309029474277013E-11,
	-0.464216300971708E-15, -0.390499637961161E-13,

	-0.236716126781431E-09, +0.454652854268717E-11,
	-0.422271787482497E-02, +0.283911742354706E-10,
	+0.270929002720228E+01
};

static const double Ipows[] = {
	0, 1, 2, 3, 4, 5, 6, 10, 12, 14, 16, 18, 20, 24, 28
};

static const int I[] = {
	0,

	0, 0, 0, 1, 1, 1, 2, 2, 3, 4, 4, 5, 5, 5, 6, 7, 8, 8, 9, 9, 9, 10,
	11, 12, 12, 13, 13, 14, 14
};

static const double Jpows[] = {
	-12, -10, -8, -6, -5, -3, -2, -1, 0, 1, 2, 3
};

static const int J[] = {
	0,

	7, 8, 9, 6, 7, 9, 7, 9, 6, 6, 10, 5, 6, 8, 11, 3, 2, 5, 1, 2, 4, 1,
	0, 0, 1, 0, 3, 0, 4
};

static const double vstar = 0.0054; /* [m³/kg] */
static const double pstar = 25; /* [MPa] */
static const double Tstar = 670; /* [K] */

double h2o_region3j_v_pT(double p, double T)
{
	double pi = p / pstar;
	double theta = T / Tstar;

	double sum = poly_value(sqrt(pi - 0.875), theta - 0.964,
			I, Ipows, 0, 15, 0,
			J, Jpows, 8, 12, 0,
			n, 29);

	return pow4(sum) * vstar;
}
