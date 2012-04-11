/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif


#include "region1.h"
#include "xmath.h"

/* Based on IF97-Rev, s. 5.2.1: The Backward Equation T(p, h) */

/* coefficient table; n[0] added for convenience */
static const double n[] = {
	+0.00000000000000E+00,

	-0.23872489924521E+03, +0.40421188637945E+03,
	+0.11349746881718E+03, -0.58457616048039E+01,
	-0.15285482413140E-03, -0.10866707695377E-05,
	-0.13391744872602E+02, +0.43211039183559E+02,
	-0.54010067170506E+02, +0.30535892203916E+02, /* [10] */

	-0.65964749423638E+01, +0.93965400878363E-02,
	+0.11573647505340E-06, -0.25858641282073E-04,
	-0.40644363084799E-08, +0.66456186191635E-07,
	+0.80670734103027E-10, -0.93477771213947E-12,
	+0.58265442020601E-14, -0.15020185953503E-16
};

static const double Ipows[] = {
	0, 1, 2, 3, 4, 5, 6
};

static const int I[] = {
	0,

	0, 0, 0, 0, 0, 0, /* [6] */
	1, 1, 1, 1, 1, 1, 1, /* [13] */
	2, 2,
	3, 3, /* [17] */
	4, 5, 6
};

static const double Jpows[] = {
	0, 1, 2, 3, 4, 6, 10, 22, 32
};

static const int J[] = {
	0,

	0, 1, 2, 5, 7, 8, /* [6] */
	0, 1, 2, 3, 4, 6, 8, /* [13] */
	6, 8,
	6, 8, /* [17] */
	8, 8, 8
};

static const double hstar = 2500; /* [kJ/kg] */

double h2o_region1_T_ph(double p, double h) /* [MPa, kJ/kg] -> [K] */
{
	double eta = h / hstar;

	return poly_value(p, eta + 1,
			I, Ipows, 0, 7, 0,
			J, Jpows, 0, 9, 0,
			n, 20);
}
