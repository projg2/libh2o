/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include "consts.h"
#include "region2.h"
#include "xmath.h"

/* Based on IF97-Rev, s. 6: Equations for Region 2 */

/* ideal-gas part coefficients */
static const double no[] = {
	+0.00000000000000E+0,

	-0.56087911283020E-2, +0.71452738081455E-1,
	-0.40710498223928E+0, +0.14240819171444E+1,
	-0.43839511319450E+1, -0.96927686500217E+1,
	+0.10086655968018E+2, -0.28408632460772E+0,
	+0.21268463753307E-1
};

/* resident part coefficients */
static const double n[] = {
	+0.00000000000000E+00,

	-0.17731742473213E-02, -0.17834862292358E-01,
	-0.45996013696365E-01, -0.57581259083432E-01,
	-0.50325278727930E-01, -0.33032641670203E-04,
	-0.18948987516315E-03, -0.39392777243355E-02,
	-0.43797295650573E-01, -0.26674547914087E-04, /* [10] */

	+0.20481737692309E-07, +0.43870667284435E-06,
	-0.32277677238570E-04, -0.15033924542148E-02,
	-0.40668253562649E-01, -0.78847309559367E-09,
	+0.12790717852285E-07, +0.48225372718507E-06,
	+0.22922076337661E-05, -0.16714766451061E-10, /* [20] */

	-0.21171472321355E-02, -0.23895741934104E+02,
	-0.59059564324270E-17, -0.12621808899101E-05,
	-0.38946842435739E-01, +0.11256211360459E-10,
	-0.82311340897998E+01, +0.19809712802088E-07,
	+0.10406965210174E-18, -0.10234747095929E-12, /* [30] */

	-0.10018179379511E-08, -0.80882908646985E-10,
	+0.10693031879409E+00, -0.33662250574171E+00,
	+0.89185845355421E-24, +0.30629316876232E-12,
	-0.42002467698208E-05, -0.59056029685639E-25,
	+0.37826947613457E-05, -0.12768608934681E-14, /* [40] */

	+0.73087610595061E-28, +0.55414715350778E-16,
	-0.94369707241210E-06
};

static const double Ipows[] = {
	0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, /* [10] */
	16, 18, 20, 21, 22, 23, 24
};

static const int I[] = {
	0,

	1, 1, 1, 1, 1,
	2, 2, 2, 2, 2, /* [10] */
	3, 3, 3, 3, 3,
	4, 4, 4, /* [18] */
	5,
	6, 6, 6, /* [22] */
	7, 7, 7,
	8, 8,
	9, /* [28] */
	10, 10, 10,
	11, 11,
	12, /* [34] */
	13, 13, 13,
	14, 15, 16, /* [40] */
	17, 17, 17
};

static const double Jpows[] = {
	0, 1, 2, 3, 4, 6, 7, 8, 10, 11, 13, /* [10] */
	14, 16, 20, 21, 25, 26, 29, 35, 36, 39, /* [20] */
	40, 48, 50, 53, 57, 58
};

static const int J[] = {
	0,

	0, 1, 2, 3, 5,
	1, 2, 4, 6, 19, /* [10] */
	0, 1, 3, 5, 18,
	1, 2, 3, /* [18] */
	6,
	3, 12, 18, /* [22] */
	0, 9, 15,
	7, 19,
	10, /* [28] */
	4, 8, 11,
	17, 23,
	25, /* [34] */
	13, 18, 22,
	14, 24, 20, /* [40] */
	16, 21, 26
};

static const double Tstar = 540; /* [K] */

static double h2o_region2_gammao_pitau(double pi, double tau, int pider, int tauder)
{
	if (!pider)
	{
		double sum = poly_value(tau, -5, 3, tauder, no);

		if (!tauder)
			sum += log(pi);

		return sum;
	}
	else
		return 1/pi;
}

static double h2o_region2_gammar_pitau(double pi, double tau, int pider, int tauder)
{
	return twoarg_poly_value(pi, tau - 0.5,
			I, Ipows, 0, 18, pider,
			J, Jpows, 0, 27, tauder,
			n, 43);
}

static double h2o_region2_gamma_pitau(double pi, double tau, int pider, int tauder)
{
	double sum;

	sum = h2o_region2_gammao_pitau(pi, tau, pider, tauder);
	sum += h2o_region2_gammar_pitau(pi, tau, pider, tauder);

	return sum;
}

double h2o_region2_v_pT(double p, double T)
{
	double pi = p;
	double tau = Tstar / T;

	double gammapi = h2o_region2_gamma_pitau(pi, tau, 1, 0);

	return pi * gammapi * R * T / p * 1E-3;
}

double h2o_region2_u_pT(double p, double T)
{
	double pi = p;
	double tau = Tstar / T;

	double gammatau = h2o_region2_gamma_pitau(pi, tau, 0, 1);
	double gammapi = h2o_region2_gamma_pitau(pi, tau, 1, 0);

	return (tau * gammatau - pi * gammapi) * R * T;
}

double h2o_region2_s_pT(double p, double T)
{
	double pi = p;
	double tau = Tstar / T;

	double gammatau = h2o_region2_gamma_pitau(pi, tau, 0, 1);
	double gamma = h2o_region2_gamma_pitau(pi, tau, 0, 0);

	return (tau * gammatau - gamma) * R;
}

double h2o_region2_h_pT(double p, double T)
{
	double pi = p;
	double tau = Tstar / T;

	double gammatau = h2o_region2_gamma_pitau(pi, tau, 0, 1);

	return tau * gammatau * R * T;
}
