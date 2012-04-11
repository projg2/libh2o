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
 * s. 4.6: Equation TB23(h,s) for Boundary between Regions 2 and 3 */

/* coefficient table; n[0] added for convenience */
static const double n[] = {
	+0.000000000000000E+00,

	+0.629096260829810E-03, -0.823453502583165E-03,
	+0.515446951519474E-07, -0.117565945784945E+01,
	+0.348519684726192E+01, -0.507837382408313E-11,
	-0.284637670005479E+01, -0.236092263939673E+01,

	+0.601492324973779E+01, +0.148039650824546E+01,
	+0.360075182221907E-03, -0.126700045009952E-01,
	-0.122184332521413E+07, +0.149276502463272E+00,
	+0.698733471798484E+00, -0.252207040114321E-01,

	+0.147151930985213E-01, -0.108618917681849E+01,
	-0.936875039816322E-03, +0.819877897570217E+02,
	-0.182041861521835E+03, +0.261907376402688E-05,
	-0.291626417025961E+05, +0.140660774926165E-04,

	+0.783237062349385E+07
};

static const double Ipows[] = {
	-12, -10, -8, -4, -3, -2, 0, 1, 3, 5, 6, 8, 12, 14
};

static const int I[] = {
	0,

	0, 1, 2, 3, 4, 5, 5, 5, 5, 6,
	7, 7, 7, 8, 8, 9, 10, 10,
	11, 11, 11, 12, 12, 13, 13
};

static const double Jpows[] = {
	-12, -8, -6, -5, -3, -2, -1, 0, 1, 2, 3, 4, 8, 10
};

static const int J[] = {
	0,

	13, 12, 10, 11, 10, 2, 9, 10, 11, 7,
	4, 5, 13, 5, 6, 3, 2, 4,
	1, 5, 6, 0, 6, 0, 8
};

static const double Tstar = 900; /* [K] */
static const double hstar = 3000; /* [kJ/kg] */
static const double sstar = 5.3; /* [kJ/kgK] */

double h2o_b23_T_hs(double h, double s)
{
	double eta = h / hstar;
	double sigma = s / sstar;

	return poly_value(eta - 0.727, sigma - 0.864,
			I, Ipows, 6, 14, 0,
			J, Jpows, 7, 14, 0,
			n, 25) * Tstar;
}
