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
 * s. 4.4: Equations h"2ab(s) and h"2c3b(s) for the Saturated Vapor Line */

/* coefficient table; n[0] added for convenience */
static const double n[] = {
	+0.00000000000000E+0,

	-0.524581170928788E03, -0.926947218142218E07,
	-0.237385107491666E03, +0.210770155812776E11,
	-0.239494562010986E02, +0.221802480294197E03,
	-0.510472533393438E07, +0.124981396109147E07,

	+0.200008436996201E10, -0.815158509791035E03,
	-0.157612685637523E03, -0.114200422332791E11,
	+0.662364680776872E16, -0.227622818296144E19,
	-0.171048081348406E32, +0.660788766938091E16,

	+0.166320055886021E23, -0.218003784381501E30,
	-0.787276140295618E30, +0.151062329700346E32,
	+0.795732170300541E07, +0.131957647355347E16,
	-0.325097068299140E24, -0.418600611419248E26,

	+0.297478906557467E35, -0.953588761745473E20,
	+0.166957699620939E25, -0.175407764869978E33,
	+0.347581490626396E35, -0.710971318427851E39
};

static const double Ipows[] = {
	0, 1, 2, 4, 7, 8, 10, 12, 18, 20, 24, 28, 32, 36
};

static const int I[] = {
	0,

	1, 1, 2, 2, 3, 3, 4, 5, 5, 6, 7, 7, 8, 9, 10,
	11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13
};

static const double Jpows[] = {
	0, 1, 2, 4, 5, 7, 8, 10, 12, 14, 20, 22, 24, 28, 32
};

static const int J[] = {
	0,

	6, 12, 3, 14, 1, 2, 5, 4, 8, 1, 0, 5, 7, 8, 14,
	6, 8, 10, 11, 12, 2, 5, 8, 9, 12, 7, 8, 10, 11, 13
};

static const double hstar = 2800; /* [kJ/kg] */
static const double s1star = 5.21; /* [kJ/kgK] */
static const double s2star = 9.2; /* [kJ/kgK] */

double h2o_b2ab4_h_s(double s)
{
	double sigma1 = s1star / s;
	double sigma2 = s / s2star;

	double sum = twoarg_poly_value(sigma1 - 0.513, sigma2 - 0.524,
			I, Ipows, 0, 14, 0,
			J, Jpows, 0, 15, 0,
			n, 30);

	return exp(sum) * hstar;
}
