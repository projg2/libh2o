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

	-0.532466612140254E+23, +0.100415480000824E+32,
	-0.191540001821367E+30, +0.105618377808847E+17,
	+0.202281884477061E+59, +0.884585472596134E+08,
	+0.166540181638363E+23, -0.313563197669111E+06,

	-0.185662327545324E+54, -0.624942093918942E-01,
	-0.504160724132590E+10, +0.187514491833092E+05,
	+0.121399979993217E-02, +0.188317043049455E+01,
	-0.167073503962060E+04, +0.965961650599775E+00,

	+0.294885696802488E+01, -0.653915627346115E+05,
	+0.604012200163444E+50, -0.198339358557937E+00,
	-0.175984090163501E+58, +0.356314881403987E+01,
	-0.575991255144384E+03, +0.456213415338071E+05,

	-0.109174044987829E+08, +0.437796099975134E+34,
	-0.616552611135792E+46, +0.193568768917797E+10,
	+0.950898170425042E+54
};

static const double Ipows[] = {
	-12, -10, -8, -6, -5, -4, -3, -2, -1, 0, 1, 3, 4, 5, 14
};

static const int I[] = {
	0,

	0, 0, 1, 2, 3, 4, 4, 5, 5, 6, 6, 7, 8, 8, 8, 9, 9, 9, 9, 10, 10, 11, 11, 11, 12, 12, 12, 13, 14
};

static const double Jpows[] = {
	0, 1, 2, 3, 4, 6, 8, 14, 16, 18, 20, 22, 24, 28, 32, 36
};

static const int J[] = {
	0,

	10, 12, 11, 7, 15, 6, 8, 5, 14, 3, 6, 4, 1, 2, 3, 0, 1, 4, 13, 0, 14, 0, 1, 2, 3, 9, 12, 4, 12
};

static const double vstar = 0.0022; /* [m³/kg] */
static const double pstar = 21; /* [MPa] */
static const double Tstar = 640; /* [K] */

double h2o_region3s_v_pT(double p, double T)
{
	double pi = p / pstar;
	double theta = T / Tstar;

	double sum = twoarg_poly_value(pi - 0.886, theta - 0.990,
			I, Ipows, 9, 15, 0,
			J, Jpows, 0, 16, 0,
			n, 29);

	return pow4(sum) * vstar;
}
