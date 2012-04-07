/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include "region4.h"
#include "xmath.h"

/* Supplementary Release on Backward Equations p(h,s) for Region 3, Equations as
 * a Function of h and s for the Region Boundaries, and an Equation Tsat(h,s)
 * for Region 4 of the IAPWS Industrial Formulation 1997 for the Thermodynamic
 * Properties of Water and Steam
 * s. 5.3: Backward Equation Tsat(h,s) */

/* coefficient table; n[0] added for convenience */
static const double n[] = {
	+0.000000000000000E+00,

	+0.179882673606601E+00, -0.267507455199603E+00,
	+0.116276722612600E+01, +0.147545428713616E+00,
	-0.512871635973248E+00, +0.421333567697984E+00,
	+0.563749522189870E+00, +0.429274443819153E+00,

	-0.335704552142140E+01, +0.108890916499278E+02,
	-0.248483390456012E+00, +0.304153221906390E+00,
	-0.494819763939905E+00, +0.107551674933261E+01,
	+0.733888415457688E-01, +0.140170545411085E-01,

	-0.106110975998808E+00, +0.168324361811875E-01,
	+0.125028363714877E+01, +0.101316840309509E+04,
	-0.151791558000712E+01, +0.524277865990866E+02,
	+0.230495545563912E+05, +0.249459806365456E-01,

	+0.210796467412137E+07, +0.366836848613065E+09,
	-0.144814105365163E+09, -0.179276373003590E-02,
	+0.489955602100459E+10, +0.471262212070518E+03,
	-0.829294390198652E+11, -0.171545662263191E+04,

	+0.355777682973575E+07, +0.586062760258436E+12,
	-0.129887635078195E+08, +0.317247449371057E+11
};

static const double Ipows[] = {
	0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 14, 16, 18, 20, 28
};

static const int I[] = {
	0,

	0, 0, 0, 1, 1, 1, 1, 2, 2, 2,
	3, 3, 3, 3, 4, 4, 5, 5, 5, 5,
	6, 6, 6, 7, 8, 8, 9, 10, 10,
	11, 11, 12, 12, 12, 13, 14
};

static const double Jpows[] = {
	0, 1, 2, 3, 4, 5, 6, 8, 12, 14, 16, 20, 22, 24, 28, 32, 36
};

static const int J[] = {
	0,

	0, 3, 8, 0, 1, 2, 5, 0, 5, 7,
	0, 2, 3, 4, 0, 1, 1, 2, 4, 10,
	6, 7, 12, 1, 11, 16, 13, 1, 14,
	8, 15, 9, 12, 16, 13, 16
};

static const double Tstar = 550; /* [K] */
static const double hstar = 2800; /* [kJ/kg] */
static const double sstar = 9.2; /* [kJ/kgK] */


double h2o_region4_T_hs(double h, double s)
{
	double eta = h / hstar;
	double etaexpr = eta - 0.119;
	double sigma = s / sstar;
	double sigmaexpr = sigma - 1.07;

	double sum = 0;

	int i;

	double etapowers[15], sigmapowers[17];

	fill_powers(etapowers, Ipows, 0, 15, etaexpr, 0);
	fill_powers(sigmapowers, Jpows, 0, 17, sigmaexpr, 0);

	for (i = 1; i <= 36; ++i)
	{
		double etapow = etapowers[I[i]];
		double sigmapow = sigmapowers[J[i]];

		sum += n[i] * etapow * sigmapow;
	}

	return sum * Tstar;
}
