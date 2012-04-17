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

	-0.401215699576099E+09, +0.484501478318406E+11,
	+0.394721471363678E-14, +0.372629967374147E+05,
	-0.369794374168666E-29, -0.380436407012452E-14,
	+0.475361629970233E-06, -0.879148916140706E-03,

	+0.844317863844331E+00, +0.122433162656600E+02,
	-0.104529634830279E+03, +0.589702771277429E+03,
	-0.291026851164444E+14, +0.170343072841850E-05,
	-0.277617606975748E-03, -0.344709605486686E+01,

	+0.221333862447095E+02, -0.194646110037079E+03,
	+0.808354639772825E-15, -0.180845209145470E-10,
	-0.696664158132412E-05, -0.181057560300994E-02,
	+0.255830298579027E+01, +0.328913873658481E+04,

	-0.173270241249904E-18, -0.661876792558034E-06,
	-0.395688923421250E-02, +0.604203299819132E-17,
	-0.400879935920517E-13, +0.160751107464958E-08,
	+0.383719409025556E-04, -0.649565446702457E-14,

	-0.149095328506000E-11, +0.541449377329581E-08
};

static const double Ipows[] = {
	-2, -1, 0, 1, 2, 5, 6, 8, 10, 12
};

static const int I[] = {
	0,

	0, 0, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4,
	4, 5, 5, 5, 6, 6, 6, 6, 7, 8, 9
};

static const double Jpows[] = {
	-12, -10, -8, -6, -5, -3, -2, -1, 0, 1, 2, 3, 4, 6, 10, 12, 14
};

static const int J[] = {
	0,

	14, 15, 4, 13, 0, 3, 6, 7, 8, 9, 10, 11, 16, 5, 6, 8, 9, 10, 2, 3,
	5, 6, 8, 12, 0, 3, 5, 0, 1, 2, 4, 0, 0, 1
};

static const double vstar = 0.0077; /* [m³/kg] */
static const double pstar = 25; /* [MPa] */
static const double Tstar = 680; /* [K] */

double h2o_region3k_v_pT(double p, double T)
{
	double pi = p / pstar;
	double theta = T / Tstar;

	return twoarg_poly_value(pi - 0.802, theta - 0.935,
			I, Ipows, 2, 10, 0,
			J, Jpows, 8, 17, 0,
			n, 34) * vstar;
}
