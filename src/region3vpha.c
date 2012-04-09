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
 * s. 3.3: Backward Equations T(p,h) and v(p,h) for Subregions 3a and 3b */

/* coefficient table; n[0] added for convenience */
static const double n[] = {
	+0.000000000000000E+0,

	+0.529944062966028E-2, -0.170099690234461E+0,
	+0.111323814312927E+2, -0.217898123145125E+4,
	-0.506061827980875E-3, +0.556495239685324E+0,
	-0.943672726094016E+1, -0.297856807561527E+0,

	+0.939353943717186E+2, +0.192944939465981E-1,
	+0.421740664704763E+0, -0.368914126282330E+7,
	-0.737566847600639E-2, -0.354753242424366E+0,
	-0.199768169338727E+1, +0.115456297059049E+1,

	+0.568366875815960E+4, +0.808169540124668E-2,
	+0.172416341519307E+0, +0.104270175292927E+1,
	-0.297691372792847E+0, +0.560394465163593E+0,
	+0.275234661176914E+0, -0.148347894866012E+0,

	-0.651142513478515E-1, -0.292468715386302E+1,
	+0.664876096952665E-1, +0.352335014263844E+1,
	-0.146340792313332E-1, -0.224503486668184E+1,
	+0.110533464706142E+1, -0.408757344495612E-1
};

static const double Ipows[] = {
	-12, -10, -8, -6, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 8
};

static const int I[] = {
	0,

	0, 0, 0, 0, 1, 1, 1, 2, 2, 3, 3, 3,
	4, 4, 5, 6, 6, 7, 7, 7, 7, 8, 8,
	9, 9, 9, 10, 10, 11, 12, 13, 14
};

static const double Jpows[] = {
	0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 16, 18, 22
};

static const int J[] = {
	0,

	6, 8, 10, 12, 4, 7, 9, 5, 10, 3, 4, 13,
	2, 3, 7, 3, 11, 0, 1, 2, 3, 0, 1,
	0, 1, 2, 0, 2, 0, 2, 2, 2
};

static const double vstar = 0.0028; /* [m³/kg] */
static const double pstar = 100; /* [MPa] */
static const double hstar = 2100; /* [kJ/kg] */

double h2o_region3a_v_ph(double p, double h)
{
	double pi = p / pstar;
	double eta = h / hstar;

	double piexpr = pi + 0.128;
	double etaexpr = eta - 0.727;

	double sum = 0;

	int i;

	double pipowers[15], etapowers[14];

	fill_powers(pipowers, Ipows, 8, 15, piexpr, 0);
	fill_powers(etapowers, Jpows, 0, 14, etaexpr, 0);

	for (i = 1; i <= 32; ++i)
	{
		double pipow = pipowers[I[i]];
		double etapow = etapowers[J[i]];

		sum += n[i] * pipow * etapow;
	}

	return sum * vstar;
}
