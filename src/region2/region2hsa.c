/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include "region2.h"
#include "xmath.h"

/* Supplementary Release on Backward Equations for Pressure as a Function
 * of Enthalpy and Entropy p(h,s) to the IAPWS Industrial Formulation 1997
 * for the Thermodynamic Properties of Water and Steam;
 * s. 6: Backward Equation p(h,s) for Region 1 */

/* coefficient table; n[0] added for convenience */
static const double n[] = {
	+0.000000000000000E+0,

	-0.182575361923032E-1, -0.125229548799536E+0,
	+0.592290437320145E+0, +0.604769706185122E+1,
	+0.238624965444474E+3, -0.298639090222922E+3,
	+0.512250813040750E-1, -0.437266515606486E+0,

	+0.413336902999504E+0, -0.516468254574773E+1,
	-0.557014838445711E+1, +0.128555037824478E+2,
	+0.114144108953290E+2, -0.119504225652714E+3,
	-0.284777985961560E+4, +0.431757846408006E+4,

	+0.112894040802650E+1, +0.197409186206319E+4,
	+0.151612444706087E+4, +0.141324451421235E-1,
	+0.585501282219601E+0, -0.297258075863012E+1,
	+0.594567314847319E+1, -0.623656565798905E+4,

	+0.965986235133332E+4, +0.681500934948134E+1,
	-0.633207286824489E+4, -0.558919224465760E+1,
	+0.400645798472063E-1
};

static const int I[] = {
	0,

	0, 0, 0, 0, 0, 0,
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
	2, 2, 2,
	3, 3, 3, 3, 3,
	4, 5, 5, 6, 7
};

static const double Jpows[] = {
	0, 1, 2, 3, 5, 6, 10, 16, 20, 22
};

static const int J[] = {
	0,

	1, 3, 5, 7, 8, 9,
	0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
	3, 7, 8,
	0, 2, 3, 5, 7,
	7, 3, 7, 3, 1
};

static const double pstar = 4; /* [MPa] */
static const double hstar = 4200; /* [kJ/kg] */
static const double sstar = 12; /* [kJ/kgK] */

double h2o_region2a_p_hs(double h, double s)
{
	double eta = h / hstar;
	double etaexpr = eta - 0.5;

	double sigma = s / sstar;
	double sigmaexpr = sigma - 1.2;

	double sum = 0;

	int i;

	double etapowers[8], sigmapowers[10];

	fill_powers_incr(etapowers, 8, etaexpr, 0);
	fill_powers(sigmapowers, Jpows, 0, 10, sigmaexpr, 0);

	for (i = 1; i <= 29; ++i)
	{
		double etapow = etapowers[I[i]];
		double sigmapow = sigmapowers[J[i]];

		sum += n[i] * etapow * sigmapow;
	}

	return pow4(sum) * pstar;
}
