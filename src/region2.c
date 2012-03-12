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
static const double no[] = {
	+0.00000000000000E+0,

	-0.96927686500217E+1,
	+0.10086655968018E+2,
	-0.56087911283020E-2,
	+0.71452738081455E-1,
	-0.40710498223928E+0,
	+0.14240819171444E+1,
	-0.43839511319450E+1,
	-0.28408632460772E+0,
	+0.21268463753307E-1
};

static const int Jo[] = {
	0,

	0, 1, -5, -4, -3, -2, -1, 2, 3
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
		for (i = 1; i <= 9; ++i)
		{
			double taupow = taupowers[Jo[i]];

			double memb = no[i] * taupow;
			if (tauder == 1)
				memb *= Jo[i];

			sum += memb;
		}

		return sum;
	}
	else
		return 1/pi;
}

static inline double h2o_region2_gamma_pT(double p, double T, int pider, int tauder)
	/* p [MPa], T [K], pider, tauder: 0/1 */
{
	double tau = Tstar / T;

	double sum;

	/* ideal-gas part */
	sum = h2o_region2_gammao_pitau(p, tau, pider, tauder);

	if (pider == 1)
		sum *= -p;
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
