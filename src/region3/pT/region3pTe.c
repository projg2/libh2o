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

	+0.715815808404721E+09, -0.114328360753449E+12,
	+0.376531002015720E-11, -0.903983668691157E-04,
	+0.665695908836252E+06, +0.535364174960127E+10,
	+0.794977402335603E+11, +0.922230563421437E+02,

	-0.142586073991215E+06, -0.111796381424162E+07,
	+0.896121629640760E+04, -0.669989239070491E+04,
	+0.451242538486834E-02, -0.339731325977713E+02,
	-0.120523111552278E+01, +0.475992667717124E+05,

	-0.266627750390341E+06, -0.153314954386524E-03,
	+0.305638404828265E+00, +0.123654999499486E+03,
	-0.104390794213011E+04, -0.157496516174308E-01,
	+0.685331118940253E+00, +0.178373462873903E+01,

	-0.544674124878910E+00, +0.204529931318843E+04,
	-0.228342359328752E+05, +0.413197481515899E+00,
	-0.341931835910405E+02
};

static const double Ipows[] = {
	-12, -10, -8, -6, -5, -4, -3, -2, -1, 0, 1, 2
};

static const int I[] = {
	0,

	0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 4, 5, 5, 6, 6, 6, 7, 7, 7, 7, 8, 9, 9, 10,
	10, 10, 11, 11
};

static const double Jpows[] = {
	0, 1, 2, 3, 4, 6, 7, 8, 10, 14, 16
};

static const int J[] = {
	0,

	9, 10, 3, 5, 8, 9, 10, 6, 7, 8, 5, 5, 2, 4, 2, 5, 6, 0, 1, 3, 4, 0, 0, 1,
	0, 4, 5, 0, 2
};

static const double vstar = 0.0032; /* [m³/kg] */
static const double pstar = 40; /* [MPa] */
static const double Tstar = 710; /* [K] */

double h2o_region3e_v_pT(double p, double T)
{
	double pi = p / pstar;
	double theta = T / Tstar;

	return poly_value(pi - 0.587, theta - 0.918,
			I, Ipows, 9, 12, 0,
			J, Jpows, 0, 11, 0,
			n, 29) * vstar;
}
