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

	+0.144165955660863E-02, -0.701438599628258E+13,
	-0.830946716459219E-16, +0.261975135368109E+00,
	+0.393097214706245E+03, -0.104334030654021E+05,
	+0.490112654154211E+09, -0.147104222772069E-03,

	+0.103602748043408E+01, +0.305308890065089E+01,
	-0.399745276971264E+07, +0.569233719593750E-11,
	-0.464923504407778E-01, -0.535400396512906E-17,
	+0.399988795693162E-12, -0.536479560201811E-06,

	+0.159536722411202E-01, +0.270303248860217E-14,
	+0.244247453858506E-07, -0.983430636716454E-05,
	+0.663513144224454E-01, -0.993456957845006E+01,
	+0.546491323528491E+03, -0.143365406393758E+05,

	+0.150764974125511E+06, -0.337209709340105E-09,
	+0.377501980025469E-08
};

static const double Ipows[] = {
	-8, -3, 0, 1, 3, 8, 10, 12, 14
};

static const int I[] = {
	0,

	0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 7, 8
};

static const double Jpows[] = {
	-12, -10, -8, -6, -5, -4, -3, -2, -1, 0, 1, 3, 4, 5, 6, 8, 14
};

static const int J[] = {
	0,

	14, 16, 6, 11, 12, 13, 15, 8, 9, 10, 13, 3, 7, 0, 1, 2, 4, 0, 1, 2, 3, 4, 5, 6, 7, 0, 0
};

static const double vstar = 0.0054; /* [m³/kg] */
static const double pstar = 23; /* [MPa] */
static const double Tstar = 650; /* [K] */

double h2o_region3r_v_pT(double p, double T)
{
	double pi = p / pstar;
	double theta = T / Tstar;

	return twoarg_poly_value(pi - 0.874, theta - 0.982,
			I, Ipows, 2, 9, 0,
			J, Jpows, 9, 17, 0,
			n, 27) * vstar;
}
