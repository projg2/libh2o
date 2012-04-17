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

	-0.820433843259950E+05, +0.473271518461586E+11,
	-0.805950021005413E-01, +0.328600025435980E+02,
	-0.356617029982490E+04, -0.172985781433335E+10,
	+0.351769232729192E+08, -0.775489259985144E+06,

	+0.710346691966018E-04, +0.993499883820274E+05,
	-0.642094171904570E+00, -0.612842816820083E+04,
	+0.232808472983776E+03, -0.142808220416837E-04,
	-0.643596060678456E-02, -0.428577227475614E+01,

	+0.225689939161918E+04, +0.100355651721510E-02,
	+0.333491455143516E+00, +0.109697576888873E+01,
	+0.961917379376452E+00, -0.838165632204598E-01,
	+0.247795908411492E+01, -0.319114969006533E+04
};

static const double Ipows[] = {
	-12, -10, -8, -6, -5, -4, -3, -2, -1, 0, 1
};

static const int I[] = {
	0,

	0, 0, 1, 1, 1, 1, 2, 3, 4, 4, 5, 5, 6, 7, 7, 7, 7, 8, 8, 8, 9, 10,
	10, 10
};

static const double Jpows[] = {
	0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12
};

static const int J[] = {
	0,

	9, 10, 6, 7, 8, 9, 8, 6, 2, 5, 3, 4, 3, 0, 1, 2, 4, 0, 1, 2, 0, 0,
	1, 3
};

static const double vstar = 0.0022; /* [m³/kg] */
static const double pstar = 23; /* [MPa] */
static const double Tstar = 650; /* [K] */

double h2o_region3q_v_pT(double p, double T)
{
	double pi = p / pstar;
	double theta = T / Tstar;

	double sum = twoarg_poly_value(pi - 0.848, theta - 0.983,
			I, Ipows, 9, 11, 0,
			J, Jpows, 0, 11, 0,
			n, 24);

	return pow4(sum) * vstar;
}
