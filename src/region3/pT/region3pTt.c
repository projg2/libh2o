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
	+0.000000000000000E00,

	+0.155287249586268E01, +0.664235115009031E01,
	-0.289366236727210E04, -0.385923202309848E13,
	-0.291002915783761E01, -0.829088246858083E12,
	+0.176814899675218E01, -0.534686695713469E09,

	+0.160464608687834E18, +0.196435366560186E06,
	+0.156637427541729E13, -0.178154560260006E01,
	-0.229746237623692E16, +0.385659001648006E08,
	+0.110554446790543E10, -0.677073830687349E14,

	-0.327910592086523E31, -0.341552040860644E51,
	-0.527251339709047E21, +0.245375640937055E24,
	-0.168776617209269E27, +0.358958955867578E29,
	-0.656475280339411E36, +0.355286045512301E39,

	+0.569021454413270E58, -0.700584546433113E48,
	-0.705772623326374E65, +0.166861176200148E53,
	-0.300475129680486E61, -0.668481295196808E51,
	+0.428432338620678E69, -0.444227367758304E72,

	-0.281396013562745E77
};

static const double Ipows[] = {
	0, 1, 2, 3, 4, 7, 10, 18, 20, 22, 24, 28, 32, 36
};

static const int I[] = {
	0,

	0, 0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 8, 9, 9, 10, 11, 12, 12, 12, 13
};

static const double Jpows[] = {
	0, 1, 3, 4, 6, 7, 8, 10, 12, 14, 16, 18, 20, 22, 24, 28, 32, 36
};

static const int J[] = {
	0,

	0, 1, 3, 8, 0, 7, 0, 4, 9, 2, 6, 0, 7, 2, 3, 5, 12, 17, 7, 8, 9, 10, 13, 11, 16, 13, 17, 14, 15, 13, 16, 17, 17
};

static const double vstar = 0.0088; /* [m³/kg] */
static const double pstar = 20; /* [MPa] */
static const double Tstar = 650; /* [K] */

double h2o_region3t_v_pT(double p, double T)
{
	double pi = p / pstar;
	double theta = T / Tstar;

	return poly_value(pi - 0.803, theta - 1.02,
			I, Ipows, 0, 14, 0,
			J, Jpows, 0, 18, 0,
			n, 33) * vstar;
}
