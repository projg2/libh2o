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
 * s. 4.5: Equation hB13(s) for Boundary between Regions 1 and 3 */

/* coefficient table; n[0] added for convenience */
static const double n[] = {
	+0.000000000000000E+00,

	+0.913965547600543E+00,
	-0.430944856041991E-04,
	+0.603235694765419E+02,
	+0.117518273082168E-17,
	+0.220000904781292E+00,
	-0.690815545851641E+02
};

static const double Ipows[] = {
	0, 1, 3, 5, 6
};

static const int I[] = {
	0,

	0, 1, 1, 2, 3, 4
};

static const double Jpows[] = {
	-12, -4, -3, -2, 0, 1, 2
};

static const int J[] = {
	0,

	4, 3, 6, 0, 1, 2
};

static const double hstar = 1700; /* [kJ/kg] */
static const double sstar = 3.8; /* [kJ/kgK] */

double h2o_b13_h_s(double s)
{
	double sigma = s / sstar;
	double sigmaexprI = sigma - 0.884;
	double sigmaexprJ = sigma - 0.864;

	double sum = 0;

	int i;

	double sigmapowersI[5], sigmapowersJ[7];

	fill_powers(sigmapowersI, Ipows, 0, 5, sigmaexprI, 0);
	fill_powers(sigmapowersJ, Jpows, 4, 7, sigmaexprJ, 0);

	for (i = 1; i <= 6; ++i)
	{
		double sigmapowI = sigmapowersI[I[i]];
		double sigmapowJ = sigmapowersJ[J[i]];

		sum += n[i] * sigmapowI * sigmapowJ;
	}

	return sum * hstar;
}
