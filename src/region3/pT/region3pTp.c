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

	-0.982825342010366E-04, +0.105145700850612E+01,
	+0.116033094095084E+03, +0.324664750281543E+04,
	-0.123592348610137E+04, -0.561403450013495E-01,
	+0.856677401640869E-07, +0.236313425393924E+03,

	+0.972503292350109E-02, -0.103001994531927E+01,
	-0.149653706199162E-08, -0.215743778861592E-04,
	-0.834452198291445E+01, +0.586602660564988E+00,
	+0.343480022104968E-25, +0.816256095947021E-05,

	+0.294985697916798E-02, +0.711730466276584E-16,
	+0.400954763806941E-09, +0.107766027032853E+02,
	-0.409449599138182E-06, -0.729121307758902E-05,
	+0.677107970938909E-08, +0.602745973022975E-07,

	-0.382323011855257E-10, +0.179946628317437E-02,
	-0.345042834640005E-03
};

static const double Ipows[] = {
	0, 1, 2, 3, 4, 6, 7, 8, 10, 12, 14, 16, 18, 20, 22, 24, 36
};

static const int I[] = {
	0,

	0, 0, 0, 0, 1, 2, 3, 3, 4, 5, 6, 6, 7, 8, 9, 9, 9, 10, 10, 10, 11, 12, 13, 14, 15, 15, 16
};

static const double Jpows[] = {
	-12, -10, -8, -6, -5, -4, -3, -2, -1, 0, 1, 2
};

static const int J[] = {
	0,

	8, 9, 10, 11, 10, 8, 6, 9, 7, 7, 4, 5, 7, 6, 0, 3, 4, 1, 2, 6, 2, 2, 1, 1, 0, 2, 0
};

static const double vstar = 0.0041; /* [m³/kg] */
static const double pstar = 23; /* [MPa] */
static const double Tstar = 650; /* [K] */

double h2o_region3p_v_pT(double p, double T)
{
	double pi = p / pstar;
	double theta = T / Tstar;

	return poly_value(sqrt(pi - 0.972), theta - 0.997,
			I, Ipows, 0, 17, 0,
			J, Jpows, 9, 12, 0,
			n, 27) * vstar;
}
