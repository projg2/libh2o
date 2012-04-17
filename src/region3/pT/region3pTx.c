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

	+0.377373741298151E+19, -0.507100883722913E+13,
	-0.103363225598860E+16, +0.184790814320773E-05,
	-0.924729378390945E-03, -0.425999562292738E+24,
	-0.462307771873973E-12, +0.107319065855767E+22,

	+0.648662492280682E+11, +0.244200600688281E+01,
	-0.851535733484258E+10, +0.169894481433592E+22,
	+0.215780222509020E-26, -0.320850551367334E+00,
	-0.382642448458610E+17, -0.275386077674421E-28,

	-0.563199253391666E+06, -0.326068646279314E+21,
	+0.397949001553184E+14, +0.100824008584757E-06,
	+0.162234569738433E+05, -0.432355225319745E+11,
	-0.592874245598610E+12, +0.133061647281106E+01,

	+0.157338197797544E+07, +0.258189614270853E+14,
	+0.262413209706358E+25, -0.920011937431142E-01,
	+0.220213765905426E-02, -0.110433759109547E+02,
	+0.847004870612087E+07, -0.592910695762536E+09,

	-0.183027173269660E-04, +0.181339603516302E+00,
	-0.119228759669889E+04, +0.430867658061468E+07
};

static const double Ipows[] = {
	-8, -6, -5, -4, -3, -1, 0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 14
};

static const int I[] = {
	0,

	0, 1, 2, 3, 3, 3, 4, 4, 5, 6, 6, 6, 7, 7, 8, 9, 9, 9, 10, 11, 11,
	11, 12, 13, 13, 13, 13, 14, 15, 15, 15, 15, 16, 16, 16, 16
};

static const double Jpows[] = {
	-12, -10, -8, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 14
};

static const int J[] = {
	0,

	19, 17, 17, 10, 11, 19, 7, 18, 14, 9, 13, 17, 1, 8, 15, 0, 9, 16,
	12, 3, 7, 10, 10, 3, 6, 10, 16, 2, 1, 2, 4, 5, 0, 1, 2, 3
};

static const double vstar = 0.0049; /* [m³/kg] */
static const double pstar = 23; /* [MPa] */
static const double Tstar = 650; /* [K] */

double h2o_region3x_v_pT(double p, double T)
{
	double pi = p / pstar;
	double theta = T / Tstar;

	return twoarg_poly_value(pi - 0.910, theta - 0.988,
			I, Ipows, 6, 17, 0,
			J, Jpows, 9, 20, 0,
			n, 36) * vstar;
}
