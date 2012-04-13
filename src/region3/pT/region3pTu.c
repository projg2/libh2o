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

	+0.122088349258355E+18, +0.104216468608488E+10,
	-0.882666931564652E+16, +0.259929510849499E+20,
	+0.222612779142211E+15, -0.878473585050085E+18,
	-0.314432577551552E+22, -0.216934916996285E+13,

	+0.159079648196849E+21, -0.339567617303423E+03,
	+0.884387651337836E+13, -0.843405926846418E+21,
	+0.114178193518022E+02, -0.122708229235641E-03,
	-0.106201671767107E+03, +0.903443213959313E+25,

	-0.693996270370852E+28, +0.648916718965575E-08,
	+0.718957567127851E+04, +0.105581745346187E-02,
	-0.651903203602581E+15, -0.160116813274676E+25,
	-0.510254294237837E-08, -0.152355388953402E+00,

	+0.677143292290144E+12, +0.276378438378930E+15,
	+0.116862983141686E-01, -0.301426947980171E+14,
	+0.169719813884840E-07, +0.104674840020929E+27,
	-0.108016904560140E+05, -0.990623601934295E-12,

	+0.536116483602738E+07, +0.226145963747881E+22,
	-0.488731565776210E-09, +0.151001548880670E-04,
	-0.227700464643920E+05, -0.781754507698846E+28
};

static const double Ipows[] = {
	-12, -10, -8, -6, -5, -3, -1, 0, 1, 2, 3, 5, 6, 8, 10, 12, 14
};

static const int I[] = {
	0,

	0, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 5, 6, 6, 6, 6, 7, 7, 8, 9, 9,
	10, 11, 11, 11, 12, 12, 13, 13, 14, 15, 15, 15, 16, 16, 16, 16
};

static const double Jpows[] = {
	-12, -10, -8, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 14
};

static const int J[] = {
	0,

	19, 17, 18, 19, 17, 18, 19, 16, 18, 13, 16, 18, 11, 8, 10, 18, 19,
	6, 10, 7, 14, 17, 4, 5, 11, 12, 4, 11, 2, 16, 5, 0, 5, 13, 0, 1, 3,
	15
};

static const double vstar = 0.0026; /* [m³/kg] */
static const double pstar = 23; /* [MPa] */
static const double Tstar = 650; /* [K] */

double h2o_region3u_v_pT(double p, double T)
{
	double pi = p / pstar;
	double theta = T / Tstar;

	return poly_value(pi - 0.902, theta - 0.988,
			I, Ipows, 7, 17, 0,
			J, Jpows, 9, 20, 0,
			n, 38) * vstar;
}
