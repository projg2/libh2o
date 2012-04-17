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

	+0.280967799943151E-38, +0.614869006573609E-30,
	+0.582238667048942E-27, +0.390628369238462E-22,
	+0.821445758255119E-20, +0.402137961842776E-14,
	+0.651718171878301E-12, -0.211773355803058E-07,

	+0.264953354380072E-02, -0.135031446451331E-31,
	-0.607246643970893E-23, -0.402352115234494E-18,
	-0.744938506925544E-16, +0.189917206526237E-12,
	+0.364975183508473E-05, +0.177274872361946E-25,

	-0.334952758812999E-18, -0.421537726098389E-08,
	-0.391048167929649E-01, +0.541276911564176E-13,
	+0.705412100773699E-11, +0.258585887897486E-08,
	-0.493111362030162E-10, -0.158649699894543E-05,

	-0.525037427886100E+00, +0.220019901729615E-02,
	-0.643064132636925E-02, +0.629154149015048E+02,
	+0.135147318617061E+03, +0.240560808321713E-06,
	-0.890763306701305E-03, -0.440209599407714E+04,

	-0.302807107747776E+03, +0.159158748314599E+04,
	+0.232534272709876E+06, -0.792681207132600E+06,
	-0.869871364662769E+11, +0.354542769185671E+12,
	+0.400849240129329E+15
};

static const double Ipows[] = {
	0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 18
};

static const int I[] = {
	0,

	0, 3, 4, 6, 7, 9, 10, 11, 12, 0, 3, 5, 6, 8, 10, 0, 3, 7, 10, 2, 3,
	4, 2, 4, 7, 4, 3, 5, 6, 0, 0, 3, 1, 0, 1, 0, 1, 0, 1
};

static const double Jpows[] = {
	-12, -10, -8, -6, -5, -4, -3, -2, -1, 0, 1, 2, 4, 5, 6
};

static const int J[] = {
	0,

	0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4,
	4, 4, 5, 6, 6, 6, 7, 8, 8, 9, 10, 10, 11, 12, 13, 14
};

static const double vstar = 0.0031; /* [m³/kg] */
static const double pstar = 23; /* [MPa] */
static const double Tstar = 650; /* [K] */

double h2o_region3n_v_pT(double p, double T)
{
	double pi = p / pstar;
	double theta = T / Tstar;

	double sum = twoarg_poly_value(pi - 0.976, theta - 0.997,
			I, Ipows, 0, 13, 0,
			J, Jpows, 9, 15, 0,
			n, 39);

	return exp(sum) * vstar;
}
