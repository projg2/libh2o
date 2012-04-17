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
 * s. 3.3: Backward Equations T(p,h) and v(p,h) for Subregions 3a and 3b */

/* coefficient table; n[0] added for convenience */
static const double n[] = {
	+0.000000000000000E+0,

	-0.225196934336318E-8, +0.140674363313486E-7,
	+0.233784085280560E-5, -0.331833715229001E-4,
	+0.107956778514318E-2, -0.271382067378863E+0,
	+0.107202262490333E+1, -0.853821329075382E+0,

	-0.215214194340526E-4, +0.769656088222730E-3,
	-0.431136580433864E-2, +0.453342167309331E+0,
	-0.507749535873652E+0, -0.100475154528389E+3,
	-0.219201924648793E+0, -0.321087965668917E+1,

	+0.607567815637771E+3, +0.557686450685932E-3,
	+0.187499040029550E+0, +0.905368030448107E-2,
	+0.285417173048685E+0, +0.329924030996098E-1,
	+0.239897419685483E+0, +0.482754995951394E+1,

	-0.118035753702231E+2, +0.169490044091791E+0,
	-0.179967222507787E-1, +0.371810116332674E-1,
	-0.536288335065096E-1, +0.160697101092520E+1
};

static const double Ipows[] = {
	-12, -8, -6, -4, -3, -2, -1, 0, 1, 2
};

static const int I[] = {
	0,

	0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2,
	3, 3, 3, 4, 4, 5, 5, 6, 6, 6, 6, 7, 8, 8, 9, 9
};

static const double Jpows[] = {
	0, 1, 2, 3, 4, 5, 6, 7, 8, 10
};

static const int J[] = {
	0,

	0, 1, 0, 1, 3, 6, 7, 8, 0, 1, 2, 5, 6, 9,
	3, 6, 9, 0, 2, 1, 2, 0, 1, 4, 5, 0, 0, 1, 2, 6
};

static const double vstar = 0.0088; /* [m³/kg] */
static const double pstar = 100; /* [MPa] */
static const double hstar = 2800; /* [kJ/kg] */

double h2o_region3b_v_ph(double p, double h)
{
	double pi = p / pstar;
	double eta = h / hstar;

	return twoarg_poly_value(pi + 0.0661, eta - 0.720,
			I, Ipows, 7, 10, 0,
			J, Jpows, 0, 10, 0,
			n, 30) * vstar;
}
