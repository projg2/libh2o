/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include "region3.h"
#include "xmath.h"

/* Revised Supplementary Release on Backward Equations for the Functions
 * T(p,h), v(p,h) and T(p,s), v(p,s) for Region 3 of the IAPWS Industrial
 * Formulation 1997 for the Thermodynamic Properties of Water and Steam
 * s. 3.4: Backward Equations T(p,s) and v(p,s) for Subregions 3a and 3b */

/* coefficient table; n[0] added for convenience */
static const double n[] = {
	+0.000000000000000E+0,

	+0.591599780322238E-4, -0.185465997137856E-2,
	+0.104190510480013E-1, +0.598647302038590E-2,
	-0.771391189901699E+0, +0.172549765557036E+1,
	-0.467076079846526E-3, +0.134533823384439E-1,

	-0.808094336805495E-1, +0.508139374365767E+0,
	+0.128584643361683E-2, -0.163899353915435E+1,
	+0.586938199318063E+1, -0.292466667918613E+1,
	-0.614076301499537E-2, +0.576199014049172E+1,

	-0.121613320606788E+2, +0.167637540957944E+1,
	-0.744135838773463E+1, +0.378168091437659E-1,
	+0.401432203027688E+1, +0.160279837479185E+2,
	+0.317848779347728E+1, -0.358362310304853E+1,

	-0.115995260446827E+7, +0.199256573577909E+0,
	-0.122270624794624E+0, -0.191449143716586E+2,
	-0.150448002905284E-1, +0.146407900162154E+2,
	-0.327477787188230E+1
};

static const double Ipows[] = {
	-12, -10, -8, -5, -4, -3, -2, 0, 1, 2
};

static const int I[] = {
	0,

	0, 0, 0, 0, 0, 0, 1, 1, 1, 1,
	2, 3, 3, 3, 4, 4, 4, 4, 5,
	6, 6, 6, 6, 6, 6, 7, 7, 7, 8, 8, 9
};

static const double Jpows[] = {
	0, 1, 2, 3, 4, 5, 6, 12
};

static const int J[] = {
	0,

	0, 1, 2, 3, 5, 6, 0, 1, 2, 4,
	0, 1, 2, 3, 0, 1, 2, 3, 1,
	0, 1, 2, 3, 4, 7, 0, 1, 2, 0, 2, 2
};

static const double vstar = 0.0088; /* [m³/kg] */
static const double pstar = 100; /* [MPa] */
static const double sstar = 5.3; /* [kJ/kgK] */

double h2o_region3b_v_ps(double p, double s)
{
	double pi = p / pstar;
	double sigma = s / sstar;

	return poly_value(pi + 0.298, sigma - 0.816,
			I, Ipows, 7, 10, 0,
			J, Jpows, 0, 8, 0,
			n, 31) * vstar;
}
