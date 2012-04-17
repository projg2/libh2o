/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include "region1.h"
#include "xmath.h"

/* Supplementary Release on Backward Equations for Pressure as a Function
 * of Enthalpy and Entropy p(h,s) to the IAPWS Industrial Formulation 1997
 * for the Thermodynamic Properties of Water and Steam;
 * s. 5: Backward Equation p(h,s) for Region 1 */

/* coefficient table; n[0] added for convenience */
static const double n[] = {
	+0.00000000000000E+0,

	-0.691997014660582E0, -0.183612548787560E2,
	-0.928332409297335E1, +0.659639569909906E2,
	-0.162060388912024E2, +0.450620017338667E3,
	+0.854680678224170E3, +0.607523214001162E4,

	+0.326487682621856E2, -0.269408844582931E2,
	-0.319947848334300E3, -0.928354307043320E3,
	+0.303634537455249E2, -0.650540422444146E2,
	-0.430991316516130E4, -0.747512324096068E3,

	+0.730000345529245E3, +0.114284032569021E4,
	-0.436407041874559E3
};

static const double Ipows[] = {
	0, 1, 2, 3, 4, 5
};

static const int I[] = {
	0,

	0, 0, 0, 0, 0, 0, 0, 0,
	1, 1, 1, 1, 2, 2, 2, 3, 4, 4, 5
};

static const double Jpows[] = {
	0, 1, 2, 4, 5, 6, 8, 10, 14
};

static const int J[] = {
	0,

	0, 1, 2, 3, 4, 5, 6, 8,
	0, 1, 3, 5, 0, 1, 7, 3, 1, 3, 0
};

static const double pstar = 100; /* [MPa] */
static const double hstar = 3400; /* [kJ/kg] */
static const double sstar = 7.6; /* [kJ/kgK] */

double h2o_region1_p_hs(double h, double s) /* [kJ/kg, kJ/kgK] -> [MPa] */
{
	double eta = h / hstar;
	double sigma = s / sstar;

	return twoarg_poly_value(eta + 0.05, sigma + 0.05,
			I, Ipows, 0, 6, 0,
			J, Jpows, 0, 9, 0,
			n, 19) * pstar;
}
