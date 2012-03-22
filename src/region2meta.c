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
static const double no_store[] = {
	-0.56087911283020E-2,
	+0.71452738081455E-1,
	-0.40710498223928E+0,
	+0.14240819171444E+1,
	-0.43839511319450E+1,
	-0.96937268393049E+1,
	+0.10087275970006E+2,
	-0.28408632460772E+0,
	+0.21268463753307E-1
};

/* shift for negative indices */
static const double *no = &no_store[5];

/* resident part coefficients */
static const double n[] = {
	+0.00000000000000E+00,

	-0.73362260186506E-2, -0.88223831943146E-1,
	-0.72334555213245E-1, -0.40813178534455E-2,
	+0.20097803380207E-2, -0.53045921898642E-1,
	-0.76190409086970E-2, -0.63498037657313E-2,

	-0.86043093028588E-1, +0.75321581522770E-2,
	-0.79238375446139E-2, -0.22888160778447E-3,
	-0.26456501482810E-2
};

static const int I[] = {
	0,

	1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 5
};

static const double Jpows[] = {
	0, 1, 2, 4, 5, 7, 9, 10, 11, 16
};

static const int J[] = {
	0,

	0, 2, 4, 8, 1, 5, 9, 3, 9, 5, 7, 6, 7
};

static const double Tstar = 540; /* [K] */

static inline double h2o_region2_meta_gammao_pitau(double pi, double tau, int pider, int tauder)
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

static inline double h2o_region2_meta_gammar_pitau(double pi, double tau, int pider, int tauder)
{
	double tauexpr = tau - 0.5;

	int i;
	double sum = 0;

	double pipowers[6];
	double taupowers[10];

	fill_powers_incr(pipowers, 6, pi, pider);
	fill_powers(taupowers, Jpows, 0, 10, tauexpr, tauder);

	for (i = 1; i <= 13; ++i)
	{
		double pipow = pipowers[I[i]];
		double taupow = taupowers[J[i]];

		double memb = n[i] * pipow * taupow;
		if (pider == 1)
			memb *= I[i];
		if (tauder == 1)
			memb *= Jpows[J[i]];

		sum += memb;
	}

	return sum;
}

static inline double h2o_region2_meta_gamma_pT(double p, double T, int pider, int tauder)
	/* p [MPa], T [K], pider, tauder: 0/1 */
{
	double tau = Tstar / T;

	double sum;

	/* ideal-gas part */
	sum = h2o_region2_meta_gammao_pitau(p, tau, pider, tauder);
	sum += h2o_region2_meta_gammar_pitau(p, tau, pider, tauder);

	if (pider == 1)
		sum *= p;
	if (tauder == 1)
		sum *= tau;

	return sum;
}

double h2o_region2_meta_v_pT(double p, double T) /* [MPa, K] -> [m³/kg] */
{
	double gammapi = h2o_region2_meta_gamma_pT(p, T, 1, 0);

	return gammapi * R * T / p * 1E-3;
}

double h2o_region2_meta_u_pT(double p, double T) /* [MPa, K] -> [kJ/kg] */
{
	double gammatau = h2o_region2_meta_gamma_pT(p, T, 0, 1);
	double gammapi = h2o_region2_meta_gamma_pT(p, T, 1, 0);

	return (gammatau - gammapi) * R * T;
}

double h2o_region2_meta_s_pT(double p, double T) /* [MPa, K] -> [kJ/kgK] */
{
	double gammatau = h2o_region2_meta_gamma_pT(p, T, 0, 1);
	double gamma = h2o_region2_meta_gamma_pT(p, T, 0, 0);

	return (gammatau - gamma) * R;
}

double h2o_region2_meta_h_pT(double p, double T) /* [MPa, K] -> [kJ/kg] */
{
	double gammatau = h2o_region2_meta_gamma_pT(p, T, 0, 1);

	return gammatau * R * T;
}
