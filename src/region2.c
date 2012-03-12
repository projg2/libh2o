/* libh2o -- steam & water properties
 * (c) 2012 Michał Górny
 * Released under the terms of the 2-clause BSD license
 */

#ifdef HAVE_CONFIG_H
#	include "config.h"
#endif

#include <math.h>

#include "consts.h"
#include "region2.h"
#include "xmath.h"

/* Based on IF97-Rev, s. 6: Equations for Region 2 */

/* ideal-gas part coefficients */
static const double no_store[] = {
	-0.56087911283020E-2,
	+0.71452738081455E-1,
	-0.40710498223928E+0,
	+0.14240819171444E+1,
	-0.43839511319450E+1,
	-0.96927686500217E+1,
	+0.10086655968018E+2,
	-0.28408632460772E+0,
	+0.21268463753307E-1
};

/* shift for negative indices */
static const double *no = &no_store[5];

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
	16, 16,
	18, /* [34] */
	20, 20, 20,
	21, 22, 23, /* [40] */
	24, 24, 24
};

static const int J[] = {
	0,

	0, 1, 2, 3, 6,
	1, 2, 4, 7, 36, /* [10] */
	0, 1, 3, 6, 35,
	1, 2, 3, /* [18] */
	7,
	3, 16, 35, /* [22] */
	0, 11, 25,
	8, 36,
	13, /* [28] */
	4, 10, 14,
	29, 50,
	57, /* [34] */
	20, 35, 48,
	21, 53, 39, /* [40] */
	26, 40, 58
};

static const double Tstar = 540; /* [K] */

static inline double h2o_region2_gammao_pitau(double pi, double tau, int pider, int tauder)
{
	if (!pider)
	{
		int i;
		double sum = 0;

		double taupowers_store[4 + 5];
			/* shift it for negative indices */
		double* taupowers = &taupowers_store[5];

		taupowers[tauder] = 1;
		taupowers[tauder + 1] = tau;

		for (i = tauder + 2; i <= 3; ++i)
			taupowers[i] = taupowers[i - 1] * tau;
		for (i = tauder - 1; i >= -5; --i)
			taupowers[i] = taupowers[i + 1] / tau;

		if (!tauder)
			sum = log(pi);

#pragma omp parallel for default(shared) private(i) reduction(+: sum)
		for (i = -5; i <= 3; ++i)
		{
			double taupow = taupowers[i];

			double memb = no[i] * taupow;
			if (tauder == 1)
				memb *= i;

			sum += memb;
		}

		return sum;
	}
	else
		return 1/pi;
}

static inline double h2o_region2_gammar_pitau(double pi, double tau, int pider, int tauder)
{
	double tauexpr = tau - 0.5;

	int i;
	double sum = 0;

#pragma omp parallel for default(shared) private(i) reduction(+: sum)
	for (i = 1; i <= 43; ++i)
	{
		double pipow = pow(pi, I[i] - pider);
		double taupow = pow(tauexpr, J[i] - tauder);

		double memb = n[i] * pipow * taupow;
		if (pider == 1)
			memb *= I[i];
		if (tauder == 1)
			memb *= J[i];

		sum += memb;
	}

	return sum;
}

static inline double h2o_region2_gamma_pT(double p, double T, int pider, int tauder)
	/* p [MPa], T [K], pider, tauder: 0/1 */
{
	double tau = Tstar / T;

	double sum;

	/* ideal-gas part */
	sum = h2o_region2_gammao_pitau(p, tau, pider, tauder);
	sum += h2o_region2_gammar_pitau(p, tau, pider, tauder);

	if (pider == 1)
		sum *= p;
	if (tauder == 1)
		sum *= tau;

	return sum;
}

double h2o_region2_v_pT(double p, double T) /* [MPa, K] -> [m³/kg] */
{
	double gammapi = h2o_region2_gamma_pT(p, T, 1, 0);

	return gammapi * R * T / p * 1E-3;
}

double h2o_region2_u_pT(double p, double T) /* [MPa, K] -> [kJ/kg] */
{
	double gammatau = h2o_region2_gamma_pT(p, T, 0, 1);
	double gammapi = h2o_region2_gamma_pT(p, T, 1, 0);

	return (gammatau - gammapi) * R * T;
}

double h2o_region2_s_pT(double p, double T) /* [MPa, K] -> [kJ/kgK] */
{
	double gammatau = h2o_region2_gamma_pT(p, T, 0, 1);
	double gamma = h2o_region2_gamma_pT(p, T, 0, 0);

	return (gammatau - gamma) * R;
}

double h2o_region2_h_pT(double p, double T) /* [MPa, K] -> [kJ/kg] */
{
	double gammatau = h2o_region2_gamma_pT(p, T, 0, 1);

	return gammatau * R * T;
}
