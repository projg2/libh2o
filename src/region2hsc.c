/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include "region2.h"
#include "xmath.h"

/* Supplementary Release on Backward Equations for Pressure as a Function
 * of Enthalpy and Entropy p(h,s) to the IAPWS Industrial Formulation 1997
 * for the Thermodynamic Properties of Water and Steam;
 * s. 6: Backward Equation p(h,s) for Region 1 */

/* coefficient table; n[0] added for convenience */
static const double n[] = {
	+0.000000000000000E+00,

	+0.112225607199012E+00, -0.339005953606712E+01,
	-0.320503911730094E+02, -0.197597305104900E+03,
	-0.407693861553446E+03, +0.132943775222331E+05,
	+0.170846839774007E+01, +0.373694198142245E+02,

	+0.358144365815434E+04, +0.423014446424664E+06,
	-0.751071025760063E+09, +0.523446127607898E+02,
	-0.228351290812417E+03, -0.960652417056937E+06,
	-0.807059292526074E+08, +0.162698017225669E+13,

	+0.772465073604171E+00, +0.463929973837746E+05,
	-0.137317885134128E+08, +0.170470392630512E+13,
	-0.251104628187308E+14, +0.317748830835520E+14,
	+0.538685623675312E+02, -0.553089094625169E+05,

	-0.102861522421405E+07, +0.204249418756234E+13,
	+0.273918446626977E+09, -0.263963146312685E+16,
	-0.107890854108088E+10, -0.296492620980124E+11,
	-0.111754907323424E+16
};

static const double Ipows[] = {
	0, 1, 2, 3, 4, 5, 6, 10, 12, 16
};

static const int I[] = {
	0,

	0, 0, 0, 0, 0, 0,
	1, 1, 1, 1, 1,
	2, 2, 2, 2, 2,
	3, 3, 3, 3, 3,
	4, 5, 5, 5, 5,
	6, 6, 7, 8, 9
};

static const double Jpows[] = {
	0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 14, 16, 18
};

static const int J[] = {
	0,

	0, 1, 2, 3, 4, 8,
	0, 2, 5, 8, 10,
	2, 3, 7, 9, 12,
	0, 5, 8, 11, 12,
	12, 1, 4, 6, 10,
	8, 12, 7, 7, 9
};

static const double pstar = 100; /* [MPa] */
static const double hstar = 3500; /* [kJ/kg] */
static const double sstar = 5.9; /* [kJ/kgK] */

double h2o_region2c_p_hs(double h, double s)
{
	double eta = h / hstar;
	double etaexpr = eta - 0.7;

	double sigma = s / sstar;
	double sigmaexpr = sigma - 1.1;

	double sum = 0;

	int i;

	double etapowers[10], sigmapowers[13];

	fill_powers(etapowers, Ipows, 0, 10, etaexpr, 0);
	fill_powers(sigmapowers, Jpows, 0, 13, sigmaexpr, 0);

	for (i = 1; i <= 31; ++i)
	{
		double etapow = etapowers[I[i]];
		double sigmapow = sigmapowers[J[i]];

		sum += n[i] * etapow * sigmapow;
	}

	return pow4(sum) * pstar;
}
