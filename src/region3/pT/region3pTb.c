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

	-0.827670470003621E-1, +0.416887126010565E+2,
	+0.483651982197059E-1, -0.291032084950276E+5,
	-0.111422582236948E+3, -0.202300083904014E-1,
	+0.294002509338515E+3, +0.140244997609658E+3,

	-0.344384158811459E+3, +0.361182452612149E+3,
	-0.140699677420738E+4, -0.202023902676481E-2,
	+0.171346792457471E+3, -0.425597804058632E+1,
	+0.691346085000334E-5, +0.151140509678925E-2,

	-0.416375290166236E-1, -0.413754957011042E+2,
	-0.506673295721637E+2, -0.572212965569023E-3,
	+0.608817368401785E+1, +0.239600660256161E+2,
	+0.122261479925384E-1, +0.216356057692938E+1,

	+0.398198903368642E+0, -0.116892827834085E+0,
	-0.102845919373532E+0, -0.492676637589284E+0,
	+0.655540456406790E-1, -0.240462535078530E+0,
	-0.269798180310075E-1, +0.128369435967012E+0
};

static const double Ipows[] = {
	-12, -10, -8, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4
};

static const int I[] = {
	0,

	0, 0, 1, 1, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 7, 8, 8, 9,
	9, 10, 10, 11, 12, 13, 13
};

static const double Jpows[] = {
	0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 14
};

static const int J[] = {
	0,

	8, 9, 7, 10, 7, 5, 6, 7, 5, 7, 8, 2, 4, 5, 0, 1, 2, 3, 5, 0, 2, 5, 0, 2, 0,
	1, 0, 2, 0, 2, 0, 1
};

static const double vstar = 0.0041; /* [m³/kg] */
static const double pstar = 100; /* [MPa] */
static const double Tstar = 860; /* [K] */

double h2o_region3b_v_pT(double p, double T)
{
	double pi = p / pstar;
	double theta = T / Tstar;

	return twoarg_poly_value(pi - 0.280, theta - 0.779,
			I, Ipows, 9, 14, 0,
			J, Jpows, 0, 11, 0,
			n, 32) * vstar;
}
