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
 * s. 3.4: Backward Equations T(p,s) and v(p,s) for Subregions 3a and 3b */

/* coefficient table; n[0] added for convenience */
static const double n[] = {
	+0.000000000000000E+00,

	+0.527111701601660E+00, -0.401317830052742E+02,
	+0.153020073134484E+03, -0.224799398218827E+04,
	-0.193993484669048E+00, -0.140467557893768E+01,
	+0.426799878114024E+02, +0.752810643416743E+00,

	+0.226657238616417E+02, -0.622873556909932E+03,
	-0.660823667935396E+00, +0.841267087271658E+00,
	-0.253717501764397E+02, +0.485708963532948E+03,
	+0.880531517490555E+03, +0.265015592794626E+07,

	-0.359287150025783E+00, -0.656991567673753E+03,
	+0.241768149185367E+01, +0.856873461222588E+00,
	+0.655143675313458E+00, -0.213535213206406E+00,
	+0.562974957606348E-02, -0.316955725450471E+15,

	-0.699997000152457E-03, +0.119845803210767E-01,
	+0.193848122022095E-04, -0.215095749182309E-04
};

static const double Ipows[] = {
	-12, -8, -6, -5, -4, -3, -2, 0, 1, 2, 3, 4, 5, 6, 8, 12, 14
};

static const int I[] = {
	0,

	0, 0, 0, 0, 1, 1, 1, 2, 2, 2,
	3, 3, 3, 3, 3, 4, 5, 5, 6, 7,
	9, 10, 11, 12, 13, 14, 15, 16
};

static const double Jpows[] = {
	0, 1, 2, 3, 4, 6, 7, 12, 24
};

static const int J[] = {
	0,

	1, 3, 4, 6, 0, 1, 3, 0, 2, 4,
	0, 1, 2, 4, 5, 7, 1, 5, 2, 0,
	1, 1, 0, 8, 0, 3, 1, 2
};

static const double Tstar = 860; /* [K] */
static const double pstar = 100; /* [MPa] */
static const double sstar = 5.3; /* [kJ/kgK] */

double h2o_region3b_T_ps(double p, double s)
{
	double pi = p / pstar;
	double sigma = s / sstar;

	double piexpr = pi + 0.760;
	double sigmaexpr = sigma - 0.818;

	double sum = 0;

	int i;

	double pipowers[17], sigmapowers[9];

	fill_powers(pipowers, Ipows, 7, 17, piexpr, 0);
	fill_powers(sigmapowers, Jpows, 0, 9, sigmaexpr, 0);

	for (i = 1; i <= 28; ++i)
	{
		double pipow = pipowers[I[i]];
		double sigmapow = sigmapowers[J[i]];

		sum += n[i] * pipow * sigmapow;
	}

	return sum * Tstar;
}
