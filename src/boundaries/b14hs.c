/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include "boundaries.h"
#include "xmath.h"

/* Supplementary Release on Backward Equations p(h,s) for Region 3, Equations as
 * a Function of h and s for the Region Boundaries, and an Equation Tsat(h,s)
 * for Region 4 of the IAPWS Industrial Formulation 1997 for the Thermodynamic
 * Properties of Water and Steam
 * s. 4.3: Equations h'1(s) and h'3a(s) for the Saturated Liquid Line */

/* coefficient table; n[0] added for convenience */
static const double n[] = {
	+0.00000000000000E+0,

	+0.332171191705237E+0, +0.611217706323496E-3,
	-0.882092478906822E+1, -0.455628192543250E+0,
	-0.263483840850452E-4, -0.223949661148062E+2,
	-0.428398660164013E+1, -0.616679338856916E+0,

	-0.146823031104040E+2, +0.284523138727299E+3,
	-0.113398503195444E+3, +0.115671380760859E+4,
	+0.395551267359325E+3, -0.154891257229285E+1,
	+0.194486637751291E+2, -0.357915139457043E+1,

	-0.335369414148819E+1, -0.664426796332460E+0,
	+0.323321885383934E+5, +0.331766744667084E+4,
	-0.223501257931087E+5, +0.573953875852936E+7,
	+0.173226193407919E+3, -0.363968822121321E-1,

	+0.834596332878346E-6, +0.503611916682674E+1,
	+0.655444787064505E+2
};

static const double Ipows[] = {
	0, 1, 2, 3, 4, 5, 7, 8, 12, 14, 16, 20, 22, 24, 28, 32
};

static const int I[] = {
	0,

	0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4,
	5, 5, 6, 7, 8, 8, 9, 9, 10,
	11, 11, 12, 13, 14, 15, 15
};


static const double Jpows[] = {
	0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 14, 16, 18, 22, 24, 28, 36
};

static const int J[] = {
	0,

	10, 16, 3, 11, 0, 5, 4, 16, 4, 11, 14,
	12, 14, 1, 4, 2, 4, 1, 13, 8,
	9, 15, 7, 3, 0, 6, 7
};

static const double hstar = 1700; /* [kJ/kg] */
static const double sstar = 3.8; /* [kJ/kgK] */

double h2o_b14_h_s(double s)
{
	double sigma = s / sstar;
	double sigmaexprI = sigma - 1.09;
	double sigmaexprJ = sigma + 0.366E-4;

	double sum = 0;

	int i;

	double sigmapowersI[16], sigmapowersJ[17];

	fill_powers(sigmapowersI, Ipows, 0, 16, sigmaexprI, 0);
	fill_powers(sigmapowersJ, Jpows, 0, 17, sigmaexprJ, 0);

	for (i = 1; i <= 27; ++i)
	{
		double sigmapowI = sigmapowersI[I[i]];
		double sigmapowJ = sigmapowersJ[J[i]];

		sum += n[i] * sigmapowI * sigmapowJ;
	}

	return sum * hstar;
}
