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

	+0.104351280732769E01, -0.227807912708513E01,
	+0.180535256723202E01, +0.420440834792042E00,
	-0.105721244834660E06, +0.436911607493884E25,
	-0.328032702839753E12, -0.678686760804270E16,

	+0.743957464645363E04, -0.356896445355761E20,
	+0.167590585186801E32, -0.355028625419105E38,
	+0.396611982166538E12, -0.414716268484468E41,
	+0.359080103867382E19, -0.116994334851995E41
};

static const double Ipows[] = {
	0, 1, 5, 6, 7, 8, 12, 16, 22, 24, 36
};

static const int I[] = {
	0,

	0, 0, 0, 1, 1, 2, 3, 4, 5, 5, 6, 7, 8, 8, 9, 10
};

static const double Jpows[] = {
	0, 1, 2, 3, 4, 7, 12, 16, 20, 32, 36
};

static const int J[] = {
	0,

	0, 3, 4, 0, 6, 10, 6, 7, 2, 8, 9, 10, 2, 9, 5, 8
};

static const double hstar = 2800; /* [kJ/kg] */
static const double sstar = 5.9; /* [kJ/kgK] */

double h2o_b2c3b4_h_s(double s)
{
	double sigma = s / sstar;

	double sum = poly_value(sigma - 1.02, sigma - 0.726,
			I, Ipows, 0, 11, 0,
			J, Jpows, 0, 11, 0,
			n, 16);

	return pow4(sum) * hstar;
}
