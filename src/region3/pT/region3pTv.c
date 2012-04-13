/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include "region3.h"
#include "xmath.h"

/* Supplementary Release on Backward Equations for Specific Volume
 * as a Function of Pressure and Temperature v(p,T)
 * for Region 3 of the IAPWS Industrial Formulation 1997 for the
 * Thermodynamic Properties of Water and Steam */

/* coefficient table; n[0] added for convenience */
static const double n[] = {
	+0.000000000000000E+00,

	-0.415652812061591E-54, +0.177441742924043E-60,
	-0.357078668203377E-54, +0.359252213604114E-25,
	-0.259123736380269E+02, +0.594619766193460E+05,
	-0.624184007103158E+11, +0.313080299915944E+17,

	+0.105006446192036E-08, -0.192824336984852E-05,
	+0.654144373749937E+06, +0.513117462865044E+13,
	-0.697595750347391E+19, -0.103977184454767E+29,
	+0.119563135540666E-47, -0.436677034051655E-41,

	+0.926990036530639E-29, +0.587793105620748E+21,
	+0.280375725094731E-17, -0.192359972440634E+23,
	+0.742705723302738E+27, -0.517429682450605E+02,
	+0.820612048645469E+07, -0.188214882341448E-08,

	+0.184587261114837E-01, -0.135830407782663E-05,
	-0.723681885626348E+17, -0.223449194054124E+27,
	-0.111526741826431E-34, +0.276032601145151E-28,
	+0.134856491567853E+15, +0.652440293345860E-09,

	+0.510655119774360E+17, -0.468138358908732E+32,
	-0.760667491183279E+16, -0.417247986986821E-18,
	+0.312545677756104E+14, -0.100375333864186E+15,
	+0.247761392329058E+27
};

static const double Ipows[] = {
	-10, -8, -6, -5, -4, -3, -2, -1, 0, 1, 3, 4, 5, 8, 10, 12, 14
};

static const int I[] = {
	0,

	0, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6,
	7, 7, 8, 8, 8, 9, 9, 10, 11, 11, 11, 12, 13, 14, 15, 16
};

static const double Jpows[] = {
	-12, -10, -8, -6, -3, -2, 0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 14
};

static const int J[] = {
	0,

	2, 0, 0, 4, 11, 12, 13, 14, 7, 8, 12, 13, 14, 16, 0, 1, 3, 14, 4,
	14, 15, 8, 10, 5, 6, 5, 12, 14, 0, 1, 9, 3, 9, 14, 8, 0, 5, 4, 7
};

static const double vstar = 0.0031; /* [m³/kg] */
static const double pstar = 23; /* [MPa] */
static const double Tstar = 650; /* [K] */

double h2o_region3v_v_pT(double p, double T)
{
	double pi = p / pstar;
	double theta = T / Tstar;

	return poly_value(pi - 0.960, theta - 0.995,
			I, Ipows, 8, 17, 0,
			J, Jpows, 6, 17, 0,
			n, 39) * vstar;
}
