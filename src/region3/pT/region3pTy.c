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

	-0.525597995024633E-09, +0.583441305228407E+04,
	-0.134778968457925E+17, +0.118973500934212E+26,
	-0.159096490904708E+27, -0.315839902302021E-06,
	+0.496212197158239E+03, +0.327777227273171E+19,

	-0.527114657850696E+22, +0.210017506281863E-16,
	+0.705106224399834E+21, -0.266713136106469E+31,
	-0.145370512554562E-07, +0.149333917053130E+28,
	-0.149795620287641E+08, -0.381881906271100E+16,

	+0.724660165585797E-04, -0.937808169550193E+14,
	+0.514411468376383E+10, -0.828198594040141E+05
};

static const double Ipows[] = {
	0, 1, 2, 3, 4, 5, 8, 10, 12
};

static const int I[] = {
	0,

	0, 0, 0, 0, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 5, 5, 6, 6, 7, 8
};

static const double Jpows[] = {
	-8, -6, -5, -4, -3, -2, -1, 0, 1, 4, 5, 6, 8
};

static const int J[] = {
	0,

	4, 8, 10, 12, 12, 3, 6, 9, 10, 0, 9, 12, 1, 11, 5, 8, 0, 5, 2, 0
};

static const double vstar = 0.0031; /* [m³/kg] */
static const double pstar = 22; /* [MPa] */
static const double Tstar = 650; /* [K] */

double h2o_region3y_v_pT(double p, double T)
{
	double pi = p / pstar;
	double theta = T / Tstar;

	double sum = poly_value(pi - 0.996, theta - 0.994,
			I, Ipows, 0, 9, 0,
			J, Jpows, 7, 13, 0,
			n, 20);

	return pow4(sum) * vstar;
}
