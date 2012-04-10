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

	+0.150042008263875E+10, -0.159397258480424E+12,
	+0.502181140217975E-03, -0.672057767855466E+02,
	+0.145058545404456E+04, -0.823889534888890E+04,
	-0.154852214233853E+00, +0.112305046746695E+02,

	-0.297000213482822E+02, +0.438565132635495E+11,
	+0.137837838635464E-02, -0.297478527157462E+01,
	+0.971777947349413E+13, -0.571527767052398E-04,
	+0.288307949778420E+05, -0.744428289262703E+14,

	+0.128017324848921E+02, -0.368275545889071E+03,
	+0.664768904779177E+16, +0.449359251958880E-01,
	-0.422897836099655E+01, -0.240614376434179E+00,
	-0.474341365254924E+01, +0.724093999126110E+00,

	+0.923874349695897E+00, +0.399043655281015E+01,
	+0.384066651868009E-01, -0.359344365571848E-02,
	-0.735196448821653E+00, +0.188367048396131E+00,
	+0.141064266818704E-03, -0.257418501496337E-02,

	+0.123220024851555E-02
};

static const double Ipows[] = {
	-12, -10, -8, -6, -5, -4, -2, -1, 0, 1, 2, 3, 8, 10
};

static const int I[] = {
	0,

	0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3,
	4, 4, 4, 5, 5, 5, 6, 6, 7, 7, 8, 8, 8,
	9, 10, 10, 11, 12, 12, 13
};

static const double Jpows[] = {
	0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 28, 32, 36
};

static const int J[] = {
	0,

	12, 13, 4, 9, 10, 11, 5, 7, 8, 12, 2, 6, 13,
	0, 11, 13, 6, 9, 14, 1, 4, 1, 6, 0, 1, 4,
	0, 0, 3, 2, 0, 1, 2
};

static const double Tstar = 760; /* [K] */
static const double pstar = 100; /* [MPa] */
static const double sstar = 4.4; /* [kJ/kgK] */

double h2o_region3a_T_ps(double p, double s)
{
	double pi = p / pstar;
	double sigma = s / sstar;

	double piexpr = pi + 0.240;
	double sigmaexpr = sigma - 0.703;

	double sum = 0;

	int i;

	double pipowers[14], sigmapowers[15];

	fill_powers(pipowers, Ipows, 8, 14, piexpr, 0);
	fill_powers(sigmapowers, Jpows, 0, 15, sigmaexpr, 0);

	for (i = 1; i <= 33; ++i)
	{
		double pipow = pipowers[I[i]];
		double sigmapow = sigmapowers[J[i]];

		sum += n[i] * pipow * sigmapow;
	}

	return sum * Tstar;
}
