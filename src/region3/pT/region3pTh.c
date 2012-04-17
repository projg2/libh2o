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

	+0.561379678887577E-01, +0.774135421587083E+10,
	+0.111482975877938E-08, -0.143987128208183E-02,
	+0.193696558764920E+04, -0.605971823585005E+09,
	+0.171951568124337E+14, -0.185461154985145E+17,

	+0.387851168078010E-16, -0.395464327846105E-13,
	-0.170875935679023E+03, -0.212010620701220E+04,
	+0.177683337348191E+08, +0.110177443629575E+02,
	-0.234396091693313E+06, -0.656174421999594E+07,

	+0.156362212977396E-04, -0.212946257021400E+01,
	+0.135249306374858E+02, +0.177189164145813E+00,
	+0.139499167345464E+04, -0.703670932036388E-02,
	-0.152011044389648E+00, +0.981916922991113E-04,

	+0.147199658618076E-02, +0.202618487025578E+02,
	+0.899345518944240E+00, -0.211346402240858E+00,
	+0.249971752957491E+02
};

static const double Ipows[] = {
	-12, -10, -8, -6, -5, -4, -3, -2, -1, 0, 1
};

static const int I[] = {
	0,

	0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 6, 6,
	7, 8, 8, 9, 10, 10
};

static const double Jpows[] = {
	0, 1, 2, 3, 4, 6, 7, 8, 10, 12, 14, 16
};

static const int J[] = {
	0,

	7, 9, 4, 5, 7, 8, 10, 11, 0, 1, 5, 6, 7, 4, 5, 7, 2, 3, 4, 2, 4, 1,
	2, 0, 0, 2, 0, 0, 2
};

static const double vstar = 0.0032; /* [m³/kg] */
static const double pstar = 25; /* [MPa] */
static const double Tstar = 660; /* [K] */

double h2o_region3h_v_pT(double p, double T)
{
	double pi = p / pstar;
	double theta = T / Tstar;

	double sum = twoarg_poly_value(pi - 0.898, theta - 0.983,
			I, Ipows, 9, 11, 0,
			J, Jpows, 0, 12, 0,
			n, 29);

	return pow4(sum) * vstar;
}
