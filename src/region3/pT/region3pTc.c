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

	+0.311967788763030E+1, +0.276713458847564E+5,
	+0.322583103403269E+8, -0.342416065095363E+3,
	-0.899732529907377E+6, -0.793892049821251E+8,
	+0.953193003217388E+2, +0.229784742345072E+4,

	+0.175336675322499E+6, +0.791214365222792E+7,
	+0.319933345844209E-4, -0.659508863555767E+2,
	-0.833426563212851E+6, +0.645734680583292E-1,
	-0.382031020570813E+7, +0.406398848470079E-4,

	+0.310327498492008E+2, -0.892996718483724E-3,
	+0.234604891591616E+3, +0.377515668966951E+4,
	+0.158646812591361E-1, +0.707906336241843E+0,
	+0.126016225146570E+2, +0.736143655772152E+0,

	+0.676544268999101E+0, -0.178100588189137E+2,
	-0.156531975531713E+0, +0.117707430048158E+2,
	+0.840143653860447E-1, -0.186442467471949E+0,
	-0.440170203949645E+2, +0.123290423502494E+7,

	-0.240650039730845E-1, -0.107077716660869E+7,
	+0.438319858566475E-1
};

static const double Ipows[] = {
	-12, -10, -8, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 8
};

static const int I[] = {
	0,

	0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 4, 4, 4, 5, 5, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 10, 10, 11, 11, 11, 11, 12, 12, 13
};

static const double Jpows[] = {
	0, 1, 2, 3, 4, 5, 6, 7, 8, 10
};

static const int J[] = {
	0,

	6, 8, 9, 6, 8, 9, 5, 6, 7, 8, 1, 4, 7, 2, 8, 0, 3, 0, 4, 5, 0, 1, 2, 0, 1, 2, 0, 2, 0, 1, 3, 7, 0, 7, 1
};

static const double vstar = 0.0022; /* [m³/kg] */
static const double pstar = 40; /* [MPa] */
static const double Tstar = 690; /* [K] */

double h2o_region3c_v_pT(double p, double T)
{
	double pi = p / pstar;
	double theta = T / Tstar;

	return twoarg_poly_value(pi - 0.259, theta - 0.903,
			I, Ipows, 9, 14, 0,
			J, Jpows, 0, 10, 0,
			n, 35) * vstar;
}
