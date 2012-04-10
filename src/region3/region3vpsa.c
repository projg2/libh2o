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
	+0.000000000000000E+00,

	+0.795544074093975E+02, -0.238261242984590E+04,
	+0.176813100617787E+05, -0.110524727080379E-02,
	-0.153213833655326E+02, +0.297544599376982E+03,
	-0.350315206871242E+08, +0.277513761062119E+00,

	-0.523964271036888E+00, -0.148011182995403E+06,
	+0.160014899374266E+07, +0.170802322663427E+13,
	+0.246866996006494E-03, +0.165326084797980E+01,
	-0.118008384666987E+00, +0.253798642355900E+01,

	+0.965127704669424E+00, -0.282172420532826E+02,
	+0.203224612353823E+00, +0.110648186063513E+01,
	+0.526127948451280E+00, +0.277000018736321E+00,
	+0.108153340501132E+01, -0.744127885357893E-01,

	+0.164094443541384E-01, -0.680468275301065E-01,
	+0.257988576101640E-01, -0.145749861944416E-03
};

static const double Ipows[] = {
	-12, -10, -8, -6, -5, -4, -3, -2, -1, 0, 1, 2, 4, 5, 6
};

static const int I[] = {
	0,

	0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2,
	3, 4, 5, 6, 6, 7, 7, 8, 8,
	9, 9, 9, 10, 11, 12, 13, 14
};

static const double Jpows[] = {
	0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 14, 16, 20, 28
};

static const int J[] = {
	0,

	8, 9, 10, 4, 7, 8, 12, 5, 6, 10, 11,
	13, 1, 5, 2, 4, 3, 7, 1, 2,
	0, 1, 3, 0, 0, 2, 2, 0
};

static const double vstar = 0.0028; /* [m³/kg] */
static const double pstar = 100; /* [MPa] */
static const double sstar = 4.4; /* [kJ/kgK] */

double h2o_region3a_v_ps(double p, double s)
{
	double pi = p / pstar;
	double sigma = s / sstar;

	double piexpr = pi + 0.187;
	double sigmaexpr = sigma - 0.755;

	double sum = 0;

	int i;

	double pipowers[15], sigmapowers[14];

	fill_powers(pipowers, Ipows, 9, 15, piexpr, 0);
	fill_powers(sigmapowers, Jpows, 0, 14, sigmaexpr, 0);

	for (i = 1; i <= 28; ++i)
	{
		double pipow = pipowers[I[i]];
		double sigmapow = sigmapowers[J[i]];

		sum += n[i] * pipow * sigmapow;
	}

	return sum * vstar;
}
