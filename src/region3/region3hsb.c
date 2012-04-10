/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include "region3.h"
#include "xmath.h"

/* Supplementary Release on Backward Equations p(h,s) for Region 3, Equations as
 * a Function of h and s for the Region Boundaries, and an Equation Tsat(h,s)
 * for Region 4 of the IAPWS Industrial Formulation 1997 for the Thermodynamic
 * Properties of Water and Steam
 * s. 3.3: Backward Equations p(h,s) */

/* coefficient table; n[0] added for convenience */
static const double n[] = {
	+0.000000000000000E+0,

	+0.125244360717979E-12, -0.126599322553713E-01,
	+0.506878030140626E+01, +0.317847171154202E+02,
	-0.391041161399932E+06, -0.975733406392044E-10,
	-0.186312419488279E+02, +0.510973543414101E+03,

	+0.373847005822362E+06, +0.299804024666572E-07,
	+0.200544393820342E+02, -0.498030487662829E-05,
	-0.102301806360030E+02, +0.552819126990325E+02,
	-0.206211367510878E+03, -0.794012232324823E+04,

	+0.782248472028153E+01, -0.586544326902468E+02,
	+0.355073647696481E+04, -0.115303107290162E-03,
	-0.175092403171802E+01, +0.257981687748160E+03,
	-0.727048374179467E+03, +0.121644822609198E-03,

	+0.393137871762692E-01, +0.704181005909296E-02,
	-0.829108200698110E+02, -0.265178818131250E+00,
	+0.137531682453991E+02, -0.522394090753046E+02,
	+0.240556298941048E+04, -0.227361631268929E+05,

	+0.890746343932567E+05, -0.239234565822486E+08,
	+0.568795808129714E+10
};

static const double Ipows[] = {
	-12, -10, -8, -6, -5, -4, -3, -2, -1, 0, 1, 2, 5, 6, 8, 10, 14
};

static const int I[] = {
	0,

	0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 3, 3, 3, 3,
	4, 5, 5, 5, 6, 6, 6, 6, 7, 7, 8, 9, 11, 11,
	12, 13, 14, 15, 16, 16
};

static const double Jpows[] = {
	0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 18, 20
};

static const int J[] = {
	0,

	2, 9, 10, 11, 13, 2, 9, 11, 12, 2, 8, 2, 6, 7, 8, 9,
	4, 5, 8, 1, 3, 5, 6, 0, 1, 0, 3, 0, 1,
	0, 1, 1, 1, 3, 7
};

static const double pstar = 16.6; /* [MPa] */
static const double hstar = 2800; /* [kJ/kg] */
static const double sstar = 5.3; /* [kJ/kgK] */

double h2o_region3b_p_hs(double h, double s)
{
	double eta = h / hstar;
	double etaexpr = eta - 0.681;

	double sigma = s / sstar;
	double sigmaexpr = sigma - 0.792;

	double sum = 0;

	int i;

	double etapowers[17], sigmapowers[14];

	fill_powers(etapowers, Ipows, 9, 17, etaexpr, 0);
	fill_powers(sigmapowers, Jpows, 0, 14, sigmaexpr, 0);

	for (i = 1; i <= 35; ++i)
	{
		double etapow = etapowers[I[i]];
		double sigmapow = sigmapowers[J[i]];

		sum += n[i] * etapow * sigmapow;
	}

	return pstar / sum;
}
