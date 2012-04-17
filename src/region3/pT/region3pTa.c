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
	+0.000000000000000E+0,

	+0.110879558823853E-2, +0.572616740810616E+3,
	-0.767051948380852E+5, -0.253321069529674E-1,
	+0.628008049345689E+4, +0.234105654131876E+6,
	+0.216867826045856E+0, -0.156237904341963E+3,

	-0.269893956176613E+5, -0.180407100085505E-3,
	+0.116732227668261E-2, +0.266987040856040E+2,
	+0.282776617243286E+5, -0.242431520029523E+4,
	+0.435217323022733E-3, -0.122494831387441E-1,

	+0.179357604019989E+1, +0.442729521058314E+2,
	-0.593223489018342E-2, +0.453186261685774E+0,
	+0.135825703129140E+1, +0.408748415856745E-1,
	+0.474686397863312E+0, +0.118646814997915E+1,

	+0.546987265727549E+0, +0.195266770452643E+0,
	-0.502268790869663E-1, -0.369645308193377E+0,
	+0.633828037528420E-2, +0.797441793901017E-1
};

static const double Ipows[] = {
	-12, -10, -8, -6, -5, -4, -3, -2, -1, 0, 1, 2
};

static const int I[] = {
	0,

	0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 4, 4, 4, 5, 6, 6, 6, 6, 7, 7, 7, 8, 8,
	8, 9, 9, 10, 10, 11, 11
};

static const double Jpows[] = {
	0, 1, 2, 3, 5, 6, 8, 10, 12
};

static const int J[] = {
	0,

	4, 7, 8, 4, 7, 8, 4, 6, 7, 1, 1, 4, 7, 6, 0, 1, 3, 5, 0, 2, 3, 0, 1,
	2, 0, 1, 0, 2, 0, 2
};

static const double vstar = 0.0024; /* [m³/kg] */
static const double pstar = 100; /* [MPa] */
static const double Tstar = 760; /* [K] */

double h2o_region3a_v_pT(double p, double T)
{
	double pi = p / pstar;
	double theta = T / Tstar;

	return twoarg_poly_value(pi - 0.085, theta - 0.817,
			I, Ipows, 9, 12, 0,
			J, Jpows, 0, 9, 0,
			n, 30) * vstar;
}
