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
	+0.000000000000000E+00,

	-0.586219133817016E-07, -0.894460355005526E+11,
	+0.531168037519774E-30, +0.109892402329239E+00,
	-0.575368389425212E-01, +0.228276853990249E+05,
	-0.158548609655002E+19, +0.329865748576503E-27,

	-0.634987981190669E-24, +0.615762068640611E-08,
	-0.961109240985747E+08, -0.406274286652625E-44,
	-0.471103725498077E-12, +0.725937724828145E+00,
	+0.187768525763682E-38, -0.103308436323771E+04,

	-0.662552816342168E-01, +0.579514041765710E+03,
	+0.237416732616644E-26, +0.271700235739893E-14,
	-0.907886213483600E+02, -0.171242509570207E-36,
	+0.156792067854621E+03, +0.923261357901470E+00,

	-0.597865988422577E+01, +0.321988767636389E+07,
	-0.399441390042203E-29, +0.493429086046981E-07,
	+0.812036983370565E-19, -0.207610284654137E-11,
	-0.340821291419719E-06, +0.542000573372233E-17,

	-0.856711586510214E-12, +0.266170454405981E-13,
	+0.858133791857099E-05
};

static const double Ipows[] = {
	-12, -10, -8, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 5, 8, 10
};

static const int I[] = {
	0,

	0, 0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 5, 5, 6, 6, 7, 7, 8, 8, 8, 9, 9,
	10, 11, 11, 12, 12, 13, 13, 13, 14, 14, 15, 15
};

static const double Jpows[] = {
	-12, -10, -8, -6, -5, -4, -3, -1, 0, 1, 2, 3, 6, 8, 14
};

static const int J[] = {
	0,

	13, 14, 7, 13, 12, 13, 14, 5, 6, 10, 13, 1, 7, 11, 1, 11, 9, 10, 2,
	5, 9, 0, 9, 7, 7, 10, 0, 4, 1, 2, 3, 0, 1, 0, 2
};

static const double vstar = 0.0039; /* [m³/kg] */
static const double pstar = 23; /* [MPa] */
static const double Tstar = 650; /* [K] */

double h2o_region3w_v_pT(double p, double T)
{
	double pi = p / pstar;
	double theta = T / Tstar;

	double sum = poly_value(pi - 0.959, theta - 0.995,
			I, Ipows, 9, 16, 0,
			J, Jpows, 8, 15, 0,
			n, 35);

	return pow4(sum) * vstar;
}
