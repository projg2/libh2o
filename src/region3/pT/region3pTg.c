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

	+0.412209020652996E-04, -0.114987238280587E+07,
	+0.948180885032080E+10, -0.195788865718971E+18,
	+0.496250704871300E+25, -0.105549884548496E+29,
	-0.758642165988278E+12, -0.922172769596101E+23,

	+0.725379072059348E+30, -0.617718249205859E+02,
	+0.107555033344858E+05, -0.379545802336487E+08,
	+0.228646846221831E+12, -0.499741093010619E+07,
	-0.280214310054101E+31, +0.104915406769586E+07,

	+0.613754229168619E+28, +0.802056715528378E+32,
	-0.298617819828065E+08, -0.910782540134681E+02,
	+0.135033227281565E+06, -0.712949383408211E+19,
	-0.104578785289542E+37, +0.304331584444093E+02,

	+0.593250797959445E+10, -0.364174062110798E+28,
	+0.921791403532461E+00, -0.337693609657471E+00,
	-0.724644143758508E+02, -0.110480239272601E+00,
	+0.536516031875059E+01, -0.291441872156205E+04,

	+0.616338176535305E+40, -0.120889175861180E+39,
	+0.818396024524612E+23, +0.940781944835829E+09,
	-0.367279669545448E+05, -0.837513931798655E+16
};

static const double Ipows[] = {
	-12, -10, -8, -6, -5, -4, -3, -2, -1, 0, 1, 3, 5, 6, 8, 10
};

static const int I[] = {
	0,

	0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4, 5, 6, 7, 7, 7, 7,
	8, 8, 8, 9, 9, 9, 10, 10, 10, 11, 12, 13, 14, 15, 15
};

static const double Jpows[] = {
	0, 1, 2, 3, 5, 6, 7, 8, 10, 12, 14, 18, 20, 22, 24
};

static const int J[] = {
	0,

	6, 9, 10, 11, 13, 14, 10, 12, 14, 6, 7, 8, 9, 7, 13, 6, 12, 13, 6,
	3, 4, 10, 14, 2, 7, 11, 0, 1, 2, 0, 1, 3, 14, 13, 9, 3, 0, 5
};


static const double vstar = 0.0027; /* [m³/kg] */
static const double pstar = 25; /* [MPa] */
static const double Tstar = 660; /* [K] */

double h2o_region3g_v_pT(double p, double T)
{
	double pi = p / pstar;
	double theta = T / Tstar;

	double sum = twoarg_poly_value(pi - 0.872, theta - 0.971,
			I, Ipows, 9, 16, 0,
			J, Jpows, 0, 15, 0,
			n, 38);

	return pow4(sum) * vstar;
}
